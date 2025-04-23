function l = tetra_edge_length( r, tetra, j )
% l = tetra_edge_length( r, tetra, j )
%
%  Length of the specified edge of the tethraedra.
%  For the edge and vertex indexing see tetra_edge_verts
%

[ i1, i2 ] = tetra_edge_verts( j );

r1 = r( tetra(:,i1), : );    
r2 = r( tetra(:,i2), : );

ev = r2 - r1;

l = sqrt( sum( ev.*ev, 2 ) );
