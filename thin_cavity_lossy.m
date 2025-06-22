function thin_cavity_lossy( show_plot )

if ~exist('show_plot', 'var')
    show_plot = 1;
end

% The geometry
w = 5e-2;  % cavity width  (x size)
l = 1e-1;  % cavity length (y size)
td = 5e-4; % dielectric thickness (plane-to-plane separation)
tc = 2e-4; % conductor thickness
to = 1e-3; % thickness of the 'outside' layers
aw = w*1.2; % simulation area width
al = l*1.2; % simulation area length
ah = td + 2*tc + 2*to; % simulation area height
x1 = -w/6;        % first port
y1 = -l/2 + l/6;
x2 =  w/4;        % second port
y2 =  l/2 - l/6;
rp =  w/80;

% Outer rectangle
foutr = @(p) drectangle(p,-aw/2,aw/2,-al/2,al/2);

% ports
fp1 = @(p) dcircle( p, x1, y1, rp );
fp2 = @(p) dcircle( p, x2, y2, rp );

% Signed distance function for (relative) mesh distribution.
% We want it denser closer to the ports.
fh = @(p) 0.001 + 0.1*abs(min(fp1(p), fp2(p)));
%% fh = @(p) ones(size(p,1),1);

h0 = w/40;

% Fixed vertices
fv = [-aw/2,-al/2;-aw/2,al/2;aw/2,-al/2;aw/2,al/2;...
      ptonrect( 0, 0, w, l, round(w/rp/4), round(l/rp/4) );...
      ptoncir( x1, y1, rp, 8 );...
      ptoncir( x2, y2, rp, 8 ); ];

[ rtri, tri ] = distmesh( foutr, fh, h0, [-aw/2,-al/2;aw/2,al/2], fv, [] );

%% patch( 'vertices', rtri, 'faces', tri, 'facecolor', [.9, .9, .9] )
%% hold on
%% scatter( fv(:,1), fv(:,2), 'r' )
%% hold off

% Number of layers
nlo = 2; % outside
nlc = 1; % conductor
nld = 1; % dielectric

zl = [ linspace( 0,  to,    nlo + 1 )(1:end-1), ...
       linspace( to, to+tc,  nlc + 1 )(1:end-1), ...
       linspace( to+tc, to+tc+td, nld + 1 )(1:end-1), ...
       linspace( to+tc+td, to+tc+td+tc, nlc + 1 )(1:end-1), ...
       linspace( to+tc+td+tc, to+tc+td+tc+to, nlo + 1 ) ];

[ r, tetra ] = extrude_mesh( rtri, tri, zl );

% Outer boundary
fout = @(p) dblock(p,-aw/2,aw/2,-al/2,al/2,0,ah);

% Conductors
fcnd1 = @(p) dblock(p,-w/2,w/2,-l/2,l/2,to,to+tc);
fcnd2 = @(p) dblock(p,-w/2,w/2,-l/2,l/2,to+tc+td,to+tc+td+tc);

for c=1:2
    % Vertices inside and on the boundary of the conductor
    fcnd = { fcnd1, fcnd2 }{ c };
    vincnd = find( fcnd(r) < 1e-7 );

    % remove tetrahedra which are inside the conductors
    tetra( ismember( tetra(:,1), vincnd ) ...
           & ismember( tetra(:,2), vincnd ) ...
           & ismember( tetra(:,3), vincnd ) ...
           & ismember( tetra(:,4), vincnd ), : ) = [ ];
end

% The mesh is now ready, generate edges.
[ edges, tetrae, tetraes ] = collect_tetra_edges( tetra );

% Surface triangles for the boundary condition
tri = surftri(r,tetra);
[ trie, tries ] = collect_tri_edges( tri, edges );

ntetra = size( tetra, 1 )
nsutri = size( tri, 1 )
nedges = size( edges, 1 )

% length of the edges
edgelen = sqrt( sum( (r(edges(:,2),:) - r(edges(:,1),:)).^2, 2) );

% Edges graph for the ports. 
G = sparse( [ edges(:,1), edges(:,2) ], ...
            [ edges(:,2), edges(:,1) ], ...
            [ edgelen,    edgelen    ], nedges, nedges );

% Find port contacts. Two ports (rows), two contacts each
port_contacts = cell( 2, 2 );

for c=1:2
    % Distance functions to find contacts vertices of the ports
    contactz = to + tc + td * (c-1);
    fport1contact = @(p) max( fp1(p), abs( p(:,3)-contactz ) );
    fport2contact = @(p) max( fp2(p), abs( p(:,3)-contactz ) );

    % Contacts vertices of the ports
    port_contacts( 1, c ) = find( abs( fport1contact( r ) ) < 1e-7 );
    port_contacts( 2, c ) = find( abs( fport2contact( r ) ) < 1e-7 );
end

% And make the ports matrix
P = mkportmat( G, edges, port_contacts );

% vertices on the outer boundary
voutb = find( abs( fout( r ) ) < 1e-7 );

% Triangles on the outer boundary, where we want to use matching bc
tout = find( ismember( tri(:,1), voutb ) ...
             & ismember( tri(:,2), voutb ) ...
             & ismember( tri(:,3), voutb ) );

% angular frequencies
freqs = linspace(1e6, 1e9, 41)*2*pi;
%% freqs = 1e9*2*pi;

Zf = [ ]; % Simulated Z for all frequency points

for w = freqs
    w

    % properties of the medium
    Z0 = sqrt(mu0/eps0);
    k0 = w*sqrt(eps0*mu0);

    % metal conductivity. Very low (copper / 100) for the loss experiments.
    sigma = 5.8e7*1e-2;

    % Full admittance matrix K
    K = sparse( nedges, nedges );

    % Tetrahedra
    for ni=1:6
        for nj=1:6
            K0 = integ_tetra_curln_curln( r, tetra, ni, nj );
            K2 = integ_tetra_n_n( r, tetra, ni, nj );
            Ke = tetraes(:,ni).*tetraes(:,nj).*( K0 - k0*k0*K2 );
            K = K + sparse( tetrae(:,ni), tetrae(:,nj), Ke, nedges, nedges );
        end
    end

    % Triangles. Surface impedance on the conductor surface and
    % matching on the outer surface.
    ktri = ones(nsutri,1)*sqrt(j*sigma*k0*Z0);
    ktri(tout) = j*k0; % matching
    for ni=1:3
        for nj=1:3
            Ks = ktri.*integ_tri_nxn_nxn( r, tri, ni, nj );
            Ke = tries(:,ni).*tries(:,nj).*Ks;
            K = K + sparse( trie(:,ni), trie(:,nj), Ke, nedges, nedges );
        end
    end

    % Rhs vector
    b = full( spdiags( -j*k0*Z0*edgelen, 0, nedges, nedges )*P );

    A = K;

    x = A \ b;

    Z = -(spdiags( edgelen, 0, nedges, nedges )*P)'*x;

    Zf = cat(3, Zf, Z);
    
end

tswrite( 'thin_cavity_highloss.z2p', freqs/(2*pi), Zf, 'Z', 50 );

if show_plot
    surft = surftri( r, tetra );

    trimesh( surft, r(:,1), r(:,2), r(:,3) );

    hold on

    %% scatter3(r(vpec,1), r(vpec,2), r(vpec,3));

    drawports( P, edges, r );

    hold off
end


