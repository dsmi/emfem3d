function [ r, tetra ] = extrude_mesh( rt, tri, zl )
% [ r, tetra ] = extrude_mesh( rt, tri, zl )
%
% 'Extrude' 2d mesh in z direction to make 3d one.
% zl is the z coordinates of the extruded layers boundaries.
%

% number of layers minus one (number of the layer iterfaces)
nl = length(zl);

% Number of vertices per layer
nvl = size( rt, 1 );

% tetrahedra mesh vertices
r = [ repmat(rt, nl, 1), ...
      reshape(repmat(reshape(zl, 1, nl), nvl, 1), nl*nvl, 1) ];

% Number of triangles per layer
ntri = size( tri, 1 );

% layer vertices offsets
layv = reshape( repmat( (0:(nl-2))*nvl, ntri, 1 ), (nl-1)*ntri, 1 );

% We need to make sure that sides of neighboring prisms have triangles
% oriented in the same way. We sort the triangle vertices by index (done
% below) and then make sure that tetrahedra of the first and second edges
% (ones which have index_v2>index_v1) are oriented one way and the third
% edge has an opposite orientation. This way tetrahedra of the triangles
% sharing and edge are always oriented in the same way.
sorted_tri = sort( tri, 2 );

v1 = repmat( sorted_tri(:,1),       nl - 1, 1 ) + layv;
v2 = repmat( sorted_tri(:,2),       nl - 1, 1 ) + layv;
v3 = repmat( sorted_tri(:,3),       nl - 1, 1 ) + layv;
v4 = repmat( sorted_tri(:,1) + nvl, nl - 1, 1 ) + layv;
v5 = repmat( sorted_tri(:,2) + nvl, nl - 1, 1 ) + layv;
v6 = repmat( sorted_tri(:,3) + nvl, nl - 1, 1 ) + layv;

tetra1 = [ v1, v2, v3, v4 ];
tetra2 = [ v2, v3, v4, v5 ];
tetra3 = [ v3, v4, v5, v6 ];
    
tetra = [ tetra1; tetra2; tetra3 ];

% Fix orientation
vol = v = tetra_v(r,tetra);
to_flip = vol<0;
tetra(to_flip,[1,2]) = tetra(to_flip,[2,1]);
