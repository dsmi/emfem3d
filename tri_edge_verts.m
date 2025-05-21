function [ i1, i2 ] = tri_edge_verts( j )
% [ i1, i2 ] = tri_edge_verts( j )
%    
%  Local vertex indices of a given edge of the triangle.
%  Edges are arranged such that edge j is 'opposite' to vertex j
%
i1e = [ 2, 3, 1 ];
i2e = [ 3, 1, 2 ];

i1 = i1e( j );
i2 = i2e( j );
