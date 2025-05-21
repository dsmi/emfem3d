function l = tri_edge_length( r, tri, j )
% l = tri_edge_length( r, tri, j )
%
%  Length of the specified edge of the triangle.
%  For the edge and vertex indexing see tri_edge_verts
%

[ i1, i2 ] = tri_edge_verts( j );

r1 = r( tri(:,i1), : );    
r2 = r( tri(:,i2), : );

ev = r2 - r1;

l = sqrt( sum( ev.*ev, 2 ) );
