
% width height and length of each bar
bw = 1e-4;
bh = 2e-4;
bl = 2e-3;

% center-to-center separation
bs = 2e-4;

% solution area dimensions
aw = 8e-4;
ah = 6e-4;
al = 3e-3;

% Outer rectangle
foutr = @(p) drectangle(p,-aw/2,aw/2,-ah/2,ah/2);

% bar rectangles
fbr1 = @(p) drectangle( p, -bs/2 - bw/2, -bs/2 + bw/2, -bh/2, bh/2 );
fbr2 = @(p) drectangle( p,  bs/2 - bw/2,  bs/2 + bw/2, -bh/2, bh/2 );

% Boundary signed distance function.
fd = @(p) foutr(p);

% Signed distance function for (relative) mesh distribution.
% We want it denser closer to the bars.
fh = @(p) 0.0002 + 0.2*min(fbr1(p), fbr2(p));

h0 = min( bw, bh ) / 3;

% fixed vertices -- corners and around each bar
fv = [ [-aw/2,-ah/2;-aw/2,ah/2;aw/2,-ah/2;aw/2,ah/2 ]; ...
       ptonrect( -bs/2, 0, bw, bh, round(bw/h0), round(bh/h0) ); ...
       ptonrect(  bs/2, 0, bw, bh, round(bw/h0), round(bh/h0) ) ];

[ rtri, tri ] = distmesh( fd, fh, h0, [-aw/2,-ah/2;aw/2,ah/2], fv, [] );

%% patch( 'vertices', rtri, 'faces', tri, 'facecolor', [.9, .9, .9] )

% length/thickness of the layer
l0 = (al-bl) / 2;

zl = [ linspace( -al/2, -bl/2, ceil( (-bl/2+al/2)/l0 ) + 1 ), ...
       linspace( -bl/2,  bl/2, ceil( bl/l0 ) + 1 )( 2:end-1 ), ...
       linspace(  bl/2,  al/2, ceil( (al/2-bl/2)/l0 ) + 1 ) ];

[ r, tetra ] = extrude_mesh( rtri, tri, zl );

% Outer boundary
fout = @(p) dblock(p,-aw/2,aw/2,-ah/2,ah/2,-al/2,al/2);

% Bars
fbar1 = @(p) dblock( p, -bs/2-bw/2, -bs/2+bw/2, -bh/2, bh/2, -bl/2, bl/2 );
fbar2 = @(p) dblock( p,  bs/2-bw/2,  bs/2+bw/2, -bh/2, bh/2, -bl/2, bl/2 );

for b=1:2
    % Vertices inside and on the boundary of the bar
    fbar = { fbar1, fbar2 }{ b };
    vinbar = find( fbar(r) < 1e-7 );

    % remove tetrahedra which are inside the bars
    tetra( ismember( tetra(:,1), vinbar ) ...
           & ismember( tetra(:,2), vinbar ) ...
           & ismember( tetra(:,3), vinbar ) ...
           & ismember( tetra(:,4), vinbar ), : ) = [ ];
end

% The mesh is now ready, generate edges.
[ edges, tetrae, tetraes ] = collect_tetra_edges( tetra );

nedges = size( edges, 1 );

% length of the edges
edgelen = sqrt( sum( (r(edges(:,2),:) - r(edges(:,1),:)).^2, 2) );

% Edges graph for the ports. 
G = sparse( [ edges(:,1), edges(:,2) ], ...
            [ edges(:,2), edges(:,1) ], ...
            [ edgelen,    edgelen    ], nedges, nedges );

% Find port contacts. Two ports (rows), two contacts each
port_contacts = cell( 2, 2 );

for m=1:2
    % Distance functions to find contacts vertices of the ports
    portz = -bl/2 + bl*(m-1);
    fc1 = @(p) dblock( p, -bs/2+bw/2, -bs/2+bw/2, -bh/2, bh/2, portz, portz );
    fc2 = @(p) dblock( p,  bs/2-bw/2,  bs/2-bw/2, -bh/2, bh/2, portz, portz );

    % Contacts vertices of the ports
    port_contacts( m, 1 ) = find( abs( fc1( r ) ) < 1e-7 );
    port_contacts( m, 2 ) = find( abs( fc2( r ) ) < 1e-7 );
end

%% % And make the ports matrix
P = mkportsmat( G, edges, port_contacts );

% vertices on the outer boundary
voutb = find( abs( fout( r ) ) < 1e-7 );

surft = surftri( r, tetra );

% remove triangles on the outer boundary
surft( ismember( surft(:,1), voutb ) ...
       & ismember( surft(:,2), voutb ) ...
       & ismember( surft(:,3), voutb ), : ) = [ ];

maxd = max( [ aw, ah, al ] );

trimesh( surft, r(:,1), r(:,2), r(:,3) );
hold on
%% scatter3(r(vc1,1), r(vc1,2), r(vc1,3), 'r')
%% scatter3(r(vc2,1), r(vc2,2), r(vc2,3), 'b')
scatter3(r(voutb,1), r(voutb,2), r(voutb,3))
%% xlim( [ -maxd maxd ] );
%% ylim( [ -maxd maxd ] );
%% zlim( [ -maxd maxd ] );

nports = size(P,2);
for j=1:nports
    [pi,pj,pw] = find(P(:,j));
    r0 = r(edges(pi,1),:);
    rr = r(edges(pi,2),:) - r(edges(pi,1),:);
    nwidx = find( pw < 0 ); % ones with neg weight are to be flipped
    r0(nwidx,:) = r(edges(pi(nwidx),2),:);
    rr(nwidx,:) = r(edges(pi(nwidx),1),:) - r(edges(pi(nwidx),2),:);
    quiver3( r0(:,1), r0(:,2), r0(:,3), rr(:,1), rr(:,2), rr(:,3), 0 );
end


hold off
