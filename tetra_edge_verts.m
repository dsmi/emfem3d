function [ i1, i2 ] = tetra_edge_verts( j )
% [ i1, i2 ] = tetra_edge_verts( j )
%    
%  Local vertex indices of a given edge of the tetrahedron.
%  Edges are arranged such that edge j is 'opposite' (does not share
%  a vertex) to edge 7-j
%
i1e = [ 1, 1, 1, 2, 4, 3 ];
i2e = [ 2, 3, 4, 3, 2, 4 ];

i1 = i1e( j );
i2 = i2e( j );
