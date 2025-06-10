
clear all;

% The geometry
l = 3e-1   % antenna length (z size)
r = 2e-3;  % radius of the antenna wire
g = r*10;  % gap between the wires at the port
al = l*2;  % simulation area length
ar = l/2;    % simulation area radius

% Outer crossection
foutx = @(p) dcircle( p, 0, 0, ar );

% wire crossection
fwire = @(p) dcircle( p, 0, 0, r );

% Signed distance function for (relative) mesh distribution.
% We want it denser closer to the wire
fh = @(p) 0.002 + 0.1*abs(fwire(p));
%% fh = @(p) ones(size(p,1),1);

h0 = r;

% Fixed vertices
fv = [ ptoncir( 0, 0, r, 8 ) ];

% 2d crossection mesh
[ xr, xtri ] = distmesh( foutx, fh, h0, [-ar/2,-ar/2;ar/2,ar/2], fv, [] );

% Number of layers
nlo = 12;  % outside
nlw = 12;  % wire (each half)
nlg = 1;   % gap

zl = [ linspace( -al/2,  -l/2,  nlo + 1 )(1:end-1), ...
       linspace(  -l/2,  -g/2,  nlw + 1 )(1:end-1), ...
       linspace(  -g/2,   g/2,  nlg + 1 )(1:end-1), ...
       linspace(   g/2,   l/2,  nlw + 1 )(1:end-1), ...
       linspace(   l/2,  al/2,  nlo + 1 ) ];

[ r, all_tetra ] = extrude_mesh( xr, xtri, zl );

% Outer boundary
fout = @(p) max( -min(-(-al/2)+p(:,3),(al/2)-p(:,3)), foutx(p) );

% wires
fwire1 = @(p) max( -min(-(-l/2)+p(:,3),(-g/2)-p(:,3)), fwire(p) );
fwire2 = @(p) max( -min(-( g/2)+p(:,3),( l/2)-p(:,3)), fwire(p) );

fwires = @(p) min( fwire1(p), fwire2(p) );

% keep tetrahedra outside the conductors
tetra = all_tetra( findtetra( r, all_tetra, @(p) -fwires(p) ), : );

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

% Find port contacts. One row for the only port.
port_contacts = cell( 1, 2 );

fcontact1 = @(p) max( fwire1(p), abs( p(:,3)+g/2 ) );
fcontact2 = @(p) max( fwire2(p), abs( p(:,3)-g/2 ) );

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
fpec = @(p) min( fwire1( p ), fwire2( p ) );

% Unknown (non-pec) edges. Check if both ends and center of the
% edge is on the pec boundary to see if this is a pec edge. 
eunk = find( ~( abs( fpec( r(edges(:,1),:) ) ) < 1e-7 ...
                & abs( fpec( r(edges(:,2),:) ) ) < 1e-7 ...
                & abs( fpec( (r(edges(:,1),:)+r(edges(:,2),:))*.5) ) < 1e-7 ) );

neunk = length(eunk);

% Matrix to get rows/columns for the unknown edges from K
S = sparse( eunk, transpose( 1:neunk ), ones( neunk, 1 ), nedges, neunk );

% Half-wavelength frequency
wh = 2*pi/(2*l*sqrt(eps0 * mu0));

% angular frequencies
%% freqs = linspace(wh/4, wh*4, 30);
freqs = wh;

Zf = [ ]; % Simulated Z for all frequency points

for w = freqs
    w

    wavelen = 2*pi/(w * sqrt(eps0 * mu0))

    % properties of the medium
    Z0 = sqrt(mu0/eps0);
    k0 = w*sqrt(eps0*mu0);

    % metal conductivity.
    sigma = 5.8e7;

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
    b = full( S'*spdiags( -j*k0*Z0*edgelen, 0, nedges, nedges )*P );

    A = S'*K*S;

    tic;
    x = A \ b;
    toc;

    Z = -(S'*spdiags( edgelen, 0, nedges, nedges )*P)'*x

    Zf = cat(3, Zf, Z);
    
end

tswrite( 'dipole_antenna.z1p', freqs/(2*pi), Zf, 'Z', 50 );

tri(tout,:) = []; % drop the outer boundary so we can see the antenna
trimesh( tri, r(:,1), r(:,2), r(:,3) );

hold on

%% %% scatter3(r(vpec,1), r(vpec,2), r(vpec,3));

drawports( P, edges, r );

patch( 'vertices', xr, 'faces', xtri, 'facecolor', [.9, .9, .9] )
%% scatter( fv(:,1), fv(:,2), 'r' )

hold off
