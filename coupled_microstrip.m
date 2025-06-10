
clear all;

% The geometry
t4 = 2e-4;   % dielectric above the trace thickness
t3 = 5e-4;   % trace thickness
t2 = 1e-3;   % trace-to-ground layer thickness
t1 = 5e-4;   % ground thickness
tw = 1e-3;   % trace width
ts = 4e-3;   % center-to-center trace separation
gw = (ts+tw)*3; % ground width
aw = gw*1.1; % simulation area width
t0 = t1;     % thickness of layer under ground
t5 = (t2+t3+t4);  % thickness of topmost layer
tl = 1e-2;   % trace length
al = tl*1.2; % simultion area length

% Crossection bounding rectangle
minx = -aw/2;
maxx =  aw/2;
miny = -t0;
maxy =  t1+t2+t3+t4+t5;

% Outer crossection
foutx = @(p) drectangle(p,minx,maxx,miny,maxy);

% traces, ground crossection and dielectric interfaces
ftrace1x = @(p) drectangle(p,-ts/2-tw/2,-ts/2+tw/2,t1+t2,t1+t2+t3);
ftrace2x = @(p) drectangle(p, ts/2-tw/2, ts/2+tw/2,t1+t2,t1+t2+t3);
fgndx = @(p) drectangle(p, -gw/2, gw/2, 0, t1 );
fditf1 = @(p) abs(-(t1+t2)+p(:,2));
fditf2 = @(p) abs(-(t1+t2+t3+t4)+p(:,2));

fh = @(p) 0.005 + 1.0*abs(min(min(min(min(ftrace1x(p), ftrace2x(p)), ...
                                      fgndx(p)), fditf1(p) ), fditf2(p)));

h0 = min( t3, tw );

linstep = @( from, to, step ) linspace(from, to, ceil((to-from)/step));
xy = @( x, y ) [ x', x'*0 + y ];

% fixed vertices -- corners, each trace, ground, dielectric interfaces
fvd = h0*0.9;
fv = [ [ minx,miny;minx,maxy;maxx,miny;maxx,maxy ]; ...
       ptonrect( -ts/2, t1+t2+t3/2, tw, t3, ceil(tw/fvd), ceil(t3/fvd) ); ...
       ptonrect(  ts/2, t1+t2+t3/2, tw, t3, ceil(tw/fvd), ceil(t3/fvd) ); ...
       ptonrect(     0,       t1/2, gw, t1, ceil(gw/fvd), ceil(t1/fvd) ); ...
       xy( linstep( -aw/2, -ts/2-tw/2, fvd )(1:end-1), t1+t2 ); ...
       xy( linstep( -ts/2+tw/2, ts/2-tw/2, fvd )(2:end-1), t1+t2 ); ...
       xy( linstep(  ts/2+tw/2, aw/2, fvd )(2:end), t1+t2 ); ...
       xy( linstep( -aw/2, aw/2, fvd )(2:end), t1+t2+t3+t4 ); ...
     ];

% 2d crossection mesh
[ xr, xtri ] = distmesh( foutx, fh, h0, [minx,miny;maxx,maxy], fv, [] );

% Number of layers
nlo = 4;   % outside at each side
nlt = 30;  % trace

zl = [ linspace( -al/2,  -tl/2,  nlo + 1 )(1:end-1), ...
       linspace( -tl/2,   tl/2,  nlt + 1 )(1:end-1), ...
       linspace(  tl/2,   al/2,  nlo + 1 ) ];

[ r, all_tetra ] = extrude_mesh( xr, xtri, zl );

% Outer boundary
fout = @(p) max( -min(-(-al/2)+p(:,3),(al/2)-p(:,3)), foutx(p) );

% traces and ground
ftrace1 = @(p) max( -min(-(-tl/2)+p(:,3),(tl/2)-p(:,3)), ftrace1x(p) );
ftrace2 = @(p) max( -min(-(-tl/2)+p(:,3),(tl/2)-p(:,3)), ftrace2x(p) );
fgnd = @(p) max( -min(-(-tl/2)+p(:,3),(tl/2)-p(:,3)), fgndx(p) );

% all conductors
fallcnd = @(p) min(min(ftrace1(p), ftrace2(p)), fgnd(p));

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
port_contacts = cell( 4, 2 );

portx = [ -ts/2, -ts/2,  ts/2, ts/2 ];
portz = [ -tl/2,  tl/2, -tl/2, tl/2 ];

for p=1:4
    px = portx(p);
    pz = portz(p);
    
    fcontact1 = @(p) dblock(p,px-tw/2,px+tw/2,t1+t2,t1+t2,pz,pz);
    fcontact2 = @(p) dblock(p,px-tw/2,px+tw/2,t1,   t1   ,pz,pz);

    % Contacts vertices of the ports
    port_contacts( p, 1 ) = find( abs( fcontact1( r ) ) < 1e-7 );
    port_contacts( p, 2 ) = find( abs( fcontact2( r ) ) < 1e-7 );
end

% And make the ports matrix
P = mkportmat( G, edges, port_contacts );

% vertices on the outer boundary
voutb = find( abs( fout( r ) ) < 1e-7 );

% Triangles on the outer boundary, where we want to use matching bc
tout = find( ismember( tri(:,1), voutb ) ...
             & ismember( tri(:,2), voutb ) ...
             & ismember( tri(:,3), voutb ) );

%% % pec boundary distance function
%% fpec = @(p) min( fwire1( p ), fwire2( p ) );

%% % Unknown (non-pec) edges. Check if both ends and center of the
%% % edge is on the pec boundary to see if this is a pec edge. 
%% eunk = find( ~( abs( fpec( r(edges(:,1),:) ) ) < 1e-7 ...
%%                 & abs( fpec( r(edges(:,2),:) ) ) < 1e-7 ...
%%                 & abs( fpec( (r(edges(:,1),:)+r(edges(:,2),:))*.5) ) < 1e-7 ) );

eunk = transpose( 1:size(edges,1) ); % no pec -- all edges unknown
neunk = length(eunk);

% Matrix to get rows/columns for the unknown edges from K
S = sparse( eunk, transpose( 1:neunk ), ones( neunk, 1 ), nedges, neunk );

% Dielectric properties for each tetrahedra
er = ones(size(tetra,1),1);
lt = zeros(size(tetra,1),1);
fref = 1e9;

% Dielectrics
fdiel1 = @(p) dblock(p,-aw/2,aw/2,-t0,   t1+t2,      -al/2,al/2);
fdiel2 = @(p) dblock(p,-aw/2,aw/2, t1+t2,t1+t2+t3+t4,-al/2,al/2);

er( findtetra( r, tetra, fdiel2 ) ) = 4.1;
lt( findtetra( r, tetra, fdiel2 ) ) = 0.026;

er( findtetra( r, tetra, fdiel1 ) ) = 3.2912;
lt( findtetra( r, tetra, fdiel1 ) ) = 0.0033;

% angular frequencies
freqs = linspace(1e6, 4e10, 40)*2*pi;
%% freqs = 1e8*2*pi;

Zf = [ ]; % Simulated Z for all frequency points

for w = freqs
    w

    % properties of the medium
    epsd = eps0*debye(er, lt, 2*pi*fref, w);
    Z0 = sqrt(mu0./epsd);
    k0 = w*sqrt(epsd*mu0);

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

    tic;
    x = A \ b;
    toc;

    Z = -(S'*spdiags( edgelen, 0, nedges, nedges )*P)'*x

    Zf = cat(3, Zf, Z);
    
end

tswrite( 'coupled_microstrip.z4p', freqs/(2*pi), Zf, 'Z', 50 );

tri(tout,:) = []; % drop the outer boundary so we can see the structure
%% trimesh( tri, r(:,1), r(:,2), r(:,3) );
trimesh( tri, r(:,3), r(:,1), r(:,2)  );

hold on

%% %% %% scatter3(r(vpec,1), r(vpec,2), r(vpec,3));

%% drawports( P, edges, r );
drawports( P, edges, [ r(:,3) r(:,1) r(:,2) ] );

%% patch( 'vertices', xr, 'faces', xtri, 'facecolor', [.9, .9, .9] )
patch('vertices', [ xr(:,1)*0 xr(:,1) xr(:,2) ], 'faces', xtri, 'facecolor', [.9 .9 .9])
%% scatter( fv(:,1), fv(:,2), 'r' )

hold off
