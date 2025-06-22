function round_cavity_matched( show_plot )

if ~exist('show_plot', 'var')
    show_plot = 1;
end

% The geometry
r1 = 5e-4;     % port radius
r2 = 1e-2;     % outer radius
t = 1e-4;      % conductor thickness
d = 5e-4;      % conductor-to-conductor separation
h = t+d+t+2*d; % full height of the simulation area

% Dielectric
er = 4.1;
lt = 0.026;
fref = 1e9;

% distance functions for the corresponding r
fr1 = @(p) sqrt((p(:,1)).^2+(p(:,2)).^2)-r1;
fr2 = @(p) sqrt((p(:,1)).^2+(p(:,2)).^2)-r2;

fh = @(p) 0.002 + 1.0*fr1(p);

h0 = r1/2;

% fixed vertices
fv = [ ptoncir( 0, 0, r1, 8 ) ];

% 2d crossection mesh
[ xr, xtri ] = distmesh( fr2, fh, h0, [-r2,-r2;r2,r2], fv, [] );

%% patch( 'vertices', xr, 'faces', xtri, 'facecolor', [.9, .9, .9] )
%% hold on
%% scatter( fv(:,1), fv(:,2), 'r' )
%% hold off

% Number of layers
nld = 1;  % between conductors
nlt = 1;  % conductor
nlo = 1;  % outside

zl = [ linspace( -h/2,     -(d+t)/2,  nlo + 1 )(1:end-1), ...
       linspace( -(d+t)/2, -d/2,      nlt + 1 )(1:end-1), ...
       linspace( -(d)/2,    d/2,      nld + 1 )(1:end-1), ...
       linspace(  d/2,     (d+t/2),   nlt + 1 )(1:end-1), ...
       linspace(  (d+t)/2,  h/2,      nlo + 1 ) ];

[ r, all_tetra ] = extrude_mesh( xr, xtri, zl );

% Outer boundary
fout = @(p) max( -min(-(-h/2)+p(:,3),(h/2)-p(:,3)), fr2(p) );

% conductors
fcnd1 = @(p) max( -min(-(-(d+t)/2)+p(:,3), (   -d/2)-p(:,3)), fr2(p) );
fcnd2 = @(p) max( -min(-(     d/2)+p(:,3), ((d+t)/2)-p(:,3)), fr2(p) );

% all conductors
fallcnd = @(p) min(fcnd1(p), fcnd2(p));

% keep tetrahedra outside the conductors
tetra = all_tetra( findtetra( r, all_tetra, @(p) -fallcnd(p) ), : );

% The mesh is now ready, generate edges.
[ edges, tetrae, tetraes ] = collect_tetra_edges( tetra );

% Surface triangles for the boundary condition
[ tri, tri_tetra ] = surftri(r,tetra);
[ trie, tries ] = collect_tri_edges( tri, edges );

ntetra = size( tetra, 1 )
nedges = size( edges, 1 )
nsurft = size( tri, 1 )

% length of the edges
edgelen = sqrt( sum( (r(edges(:,2),:) - r(edges(:,1),:)).^2, 2) );

% Edges graph for the ports. 
G = sparse( [ edges(:,1), edges(:,2) ], ...
            [ edges(:,2), edges(:,1) ], ...
            [ edgelen,    edgelen    ], nedges, nedges );

% Find port contacts.
port_contacts = cell( 1, 2 );
    
fcontact1 = @(p) max( fr1(p), abs( p(:,3)-(-d/2) ) );
fcontact2 = @(p) max( fr1(p), abs( p(:,3)-( d/2) ) );

% Contacts vertices of the ports
port_contacts( 1, 1 ) = find( abs( fcontact1( r ) ) < 1e-7 );
port_contacts( 1, 2 ) = find( abs( fcontact2( r ) ) < 1e-7 );

% And make the ports matrix
P = mkportmat( G, edges, port_contacts );

% vertices on the outer boundary
voutb = find( abs( fout( r ) ) < 1e-7 );

% Triangles on the outer boundary, where we want to use matching bc
tout = find( ismember( tri(:,1), voutb ) ...
             & ismember( tri(:,2), voutb ) ...
             & ismember( tri(:,3), voutb ) );

% pec boundary distance function
fpec = @(p) fallcnd(p);

%% % Unknown (non-pec) edges. Check if both ends and center of the
%% % edge is on the pec boundary to see if this is a pec edge. 
%% eunk = find( ~( abs( fpec( r(edges(:,1),:) ) ) < 1e-10 ...
%%                 & abs( fpec( r(edges(:,2),:) ) ) < 1e-10 ...
%%                 & abs( fpec( (r(edges(:,1),:)+r(edges(:,2),:))*.5) ) < 1e-10 ) );

eunk = transpose( 1:size(edges,1) ); % no pec -- all edges unknown

neunk = length(eunk);

% Matrix to get rows/columns for the unknown edges from K
S = sparse( eunk, transpose( 1:neunk ), ones( neunk, 1 ), nedges, neunk );

% angular frequencies
freqs = linspace(1e7, 1.2e10, 21)*2*pi;

Zf = [ ];  % Simulated Z for all frequency points
Zaf = [ ]; % Analytical Z for all frequency points
Zmf = [ ]; % Matched cavity Z for all frequency points

for w = freqs
    
    freq_hz = (w/(2*pi))

    % properties of the medium
    epsd = eps0*debye(er, lt, 2*pi*fref, w);
    Z0 = sqrt(mu0./epsd)*ones(size(tetra,1),1);
    k0 = w*sqrt(epsd*mu0)*ones(size(tetra,1),1);

    % the same over the entire domain
    k0_Z0 = w*mu0;

    % metal conductivity.
    sigma = 5.8e7;

    % Full admittance matrix K
    K = sparse( nedges, nedges );

    % Tetrahedra
    for ni=1:6
        for nj=1:6
            K0 = integ_tetra_curln_curln( r, tetra, ni, nj );
            K2 = integ_tetra_n_n( r, tetra, ni, nj );
            Ke = tetraes(:,ni).*tetraes(:,nj).*( K0 - k0.*k0.*K2 );
            K = K + sparse( tetrae(:,ni), tetrae(:,nj), Ke, nedges, nedges );
        end
    end

    % Triangles. Surface impedance on the conductor surface and
    % matching on the outer surface.
    ktri = ones(size(tri,1),1)*sqrt(j*sigma*k0_Z0);
    ktri(tout) = j*k0(tri_tetra(tout)); % matching
    for ni=1:3
        for nj=1:3
            Ks = ktri.*integ_tri_nxn_nxn( r, tri, ni, nj );
            Ke = tries(:,ni).*tries(:,nj).*Ks;
            K = K + sparse( trie(:,ni), trie(:,nj), Ke, nedges, nedges );
        end
    end

    % Rhs vector
    b = full( S'*spdiags( -j*k0_Z0*edgelen, 0, nedges, nedges )*P );

    A = S'*K*S;

    x = A \ b;

    Z = -(S'*spdiags( edgelen, 0, nedges, nedges )*P)'*x;

    Zf = cat(3, Zf, Z);

    % Analytical solution
    Yplane = j*w*epsd/d;
    Zs = 2*sqrt(j*w*mu0/sigma); % surface impedance
    Zplane = Zs+j*w*mu0*d;
    k = sqrt(-Yplane*Zplane);

    wavelen = 2*pi/k

    % besselh(0,2) is an outward wave from center, we calculate
    % current and voltage at r1
    V = besselh(0,2,k*r1);
    I = -(-besselh(1,2,k*r1)*k*r1)*2*pi/Zplane;
    Za = V(1)/I;
    Zaf = cat(3, Zaf, Za);

    tswrite( 'round_cavity_matched.z1p', freqs/(2*pi),   Zf,  'Z', 50 );
    tswrite( 'round_cavity_matched_analytical.z1p', freqs/(2*pi), Zaf, 'Z', 50 );
    
end

if show_plot
    tri(tout,:) = []; % drop the outer boundary so we can see the structure
    trimesh( tri, r(:,1), r(:,2), r(:,3) );

    hold on

    %% scatter3(r(vpec,1), r(vpec,2), r(vpec,3));

    drawports( P, edges, r );

    patch( 'vertices', xr, 'faces', xtri, 'facecolor', [.9, .9, .9] )
    scatter3( fv(:,1), fv(:,2), fv(:,2)*0, 'r' )

    hold off
end

