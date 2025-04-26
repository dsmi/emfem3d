function test_extrude_mesh()

fd = @(p) drectangle(p,-0.2,0.2,-0.3,0.3);
fh = @(p) ones(size(p,1),1);

box = [ -0.2,-0.3;0.2,0.3 ];
fix = [ -0.2,-0.3;-0.2,0.3;0.2,-0.3;0.2,0.3 ];

[rt,tri] = distmesh(fd,fh,0.2,box,fix);

% rt = [ 0 0 ; 0 1 ; 1 0 ];
% tri = [ 1 2 3 ];

% trimesh( tri, r(:,1), r(:,2)  )

[ r, tetra ] = extrude_mesh( rt, tri, [ 0, 0.2, 0.5 ] );

% tetramesh(tetra,r);

vol = tetra_v(r,tetra);

assert( min( min( vol, 1 ), 2 ) > 5e-4 )


