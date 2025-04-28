
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
fh = @(p) 0.0002 + 0.3*min(fbr1(p), fbr2(p));

h0 = min( bw, bh );

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

% vertices on the outer boundary
voutb = find( abs( fout( r ) ) < 1e-7 );

surft = surftri( r, tetra );

% remove triangles on the outer boundary
surft( ismember( surft(:,1), voutb ) ...
       & ismember( surft(:,2), voutb ) ...
       & ismember( surft(:,3), voutb ), : ) = [ ];

maxd = max( [ aw, ah, al ] );

trimesh( surft, r(:,1), r(:,2), r(:,3) );
xlim( [ -maxd maxd ] );
ylim( [ -maxd maxd ] );
zlim( [ -maxd maxd ] );
