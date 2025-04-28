function test_extrude_mesh()

fd = @(p) drectangle(p,-0.2,0.2,-0.3,0.3);
fh = @(p) ones(size(p,1),1);

box = [ -0.2,-0.3;0.2,0.3 ];
fix = [ -0.2,-0.3;-0.2,0.3;0.2,-0.3;0.2,0.3 ];

[rt,tri] = distmesh(fd,fh,0.1,box,fix);

% rt = [ 0 0 ; 0 1 ; 1 0 ];
% tri = [ 1 2 3 ];

% trimesh( tri, r(:,1), r(:,2)  )

[ r, tetra ] = extrude_mesh( rt, tri, [ 0, 0.2, 0.5 ] );

% tetramesh(tetra,r);

vol = tetra_v(r,tetra);

assert( min( vol ) > 1e-4 )

% Next, we want to make sure that internal faces of the neighboring
% tetrahedra coincide. We call surftri which removes the inner triangles
% by checking for duplicates and keeping only the unique ones, and make
% sure that only the surface ones remain.

stri = surftri(r,tetra);

fdbou = @(p) dblock(p,-0.2,0.2,-0.3,0.3,0,0.5);

% Indices of the boundary vertices
bv = find( abs( fdbou(r) ) < 1e-7 );

% all stri should be connected between the boundary vertices
btri = find( ismember( stri(:,1), bv ) ...
             & ismember( stri(:,2), bv ) ...
             & ismember( stri(:,3), bv ) );

assert( btri' == 1:size(stri,1) );

