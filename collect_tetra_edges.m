function [ edges, tetrae, tetraes ] = collect_tetra_edges(tetra)
% [ edges, tetrae, tetraes ] = collect_tetra_edges(tetra)
%
% Builds edges for a given set of tegrahedra. Each edge joins two vertices,
% each tetrahedron has six edges. (see tetra_edge_verts vertices that each
% edge connects)
% The function returns array of edges and indices and signs of the edges
% for each terahedron.
% 
%  Input:
%    tetra   - ntetra-by-4 array, vertices of the tetrahedra
%  Outputs:
%    edges   - edges, num_of_edges-by-2 array, two vertex indices
%              for each edge. edges(0,:) < edges(1,:) and edges are
%              sorted.
%    tetrae  - indices of edges of of each tetrahedron, ntetra-by-6.
%    tetraes - signs of the edges in tetrae, ntetra-by-6 array.
%              To remove duplicates we change orientation of some edges
%              (swap vertices), the sign is 1 if the tetrahedron edge
%              has the same direction as the corresponding edge in the
%              edges array, and -1 if the direction is opposite.
%

ntetra = size(tetra, 1);

% All edges - ones different by direction only not merged yet.
all_edges = zeros(ntetra*6, 2);

% Indices of the tetrahedron edges in all_edges array.
tetra_all_edges = zeros(ntetra, 6);

for j=1:6
    [ i1, i2 ] = tetra_edge_verts( j );
    tetra_all_edges(:, j) = j:6:(ntetra*6);
    all_edges(tetra_all_edges(:, j), :) = [ tetra(:,i1), tetra(:,i2) ];
end

% Unify the edges direction - vertex with the lesser index comes first.
to_flip = find(all_edges(:,1) > all_edges(:,2));
flipev1 = all_edges(to_flip,1);
all_edges(to_flip,1) = all_edges(to_flip,2);
all_edges(to_flip,2) = flipev1;

% Signs of the edges in all_edges ( 1 or -1 )
all_signs = ones(size(all_edges, 1), 1);
all_signs(to_flip) = -1;

% edges == all_edges(a2e,:); all_edges == edges(e2a,:)
[ edges, a2e, e2a ] = unique( all_edges, 'rows' );

% Now fill tetrahedron edges and signs
tetrae = zeros( ntetra, 6 );
tetraes = zeros( ntetra, 6 );

for j=1:6
    tetrae (:,j) =       e2a( tetra_all_edges(:, j) );
    tetraes(:,j) = all_signs( tetra_all_edges(:, j) );
end
