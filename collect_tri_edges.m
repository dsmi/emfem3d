function [ trie, tries ] = collect_tri_edges( tri, edges )
% [ trie, tries ] = collect_tri_edges( tri, edges )
%
% Identifies subset of edges which belong to the given set of triangles,
% and determine sign of each edge: 1 or the edge direction coincides
% with the direction of the triangle edge and -1 otherwise.
% We use this function to determine edges which belong to surface triangles
% of the tetrahedral mesh.
% 
%  Input:
%    tri     - ntri-by-3 array, vertices of the triangles
%    edges   - set of edges, which must include all the triangle
%              edges. Vetices in each edge must be sorted ascending
%              (edge(1,k)<(edge(2,k)) and therefore edge could be oriented
%              differently than the triangle edge, in this case tries
%              will be -1.
%  Outputs:
%    trie    - indices of edges of each triangle, ntri-by-3.
%    tries   - signs of the edges in trie, ntri-by-3 array.
%              If the edge in the 'edges' array is oriented the
%              same as the triangle edge it is 1, and -1 otherwise.
%

ntri = size(tri, 1);

% All edges - ones different by direction only not merged yet.
tri_edges = zeros(ntri*3, 2);

% Indices of the triangle edges in tri_edges array.
tri_edges_idx = zeros(ntri, 3);

for j=1:3
    [ i1, i2 ] = tri_edge_verts( j );
    tri_edges_idx(:, j) = j:3:(ntri*3);
    tri_edges(tri_edges_idx(:, j), :) = [ tri(:,i1), tri(:,i2) ];
end

% Unify the edges direction - vertex with the lesser index comes first.
to_flip = find(tri_edges(:,1) > tri_edges(:,2));
flipev1 = tri_edges(to_flip,1);
tri_edges(to_flip,1) = tri_edges(to_flip,2);
tri_edges(to_flip,2) = flipev1;

% Signs of the edges in tri_edges ( 1 or -1 )
tri_edges_signs = ones(size(tri_edges, 1), 1);
tri_edges_signs(to_flip) = -1;

[ edge_found, edge_idx ] = ismember(tri_edges, edges, "rows");

assert( edge_found );

% Now fill triangle edges and signs
trie = zeros( ntri, 3 );
tries = zeros( ntri, 3 );

for j=1:3
    trie (:,j) =        edge_idx( tri_edges_idx(:, j) );
    tries(:,j) = tri_edges_signs( tri_edges_idx(:, j) );
end
