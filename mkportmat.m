function P = mkportmat( G, edges, port_contacts )
% P = mkportmat( G, edges, port_contacts )
%
% Makes nedges-by-nports sparse matrix P where element p(i,j) is weight of
% edge i if it belongs to port j and 0 otherwise. Weight of the edge is
% 1/num_port_paths where num_port_paths depends on the number of contact
% vertices. Weight can be positive or negative, depending on whether
% the direction of the edge coincides with the direction of the port path.
%
% G is the adjacency matrix of the edges graph, and port_contacts is
% nports-by-2 cell array with the contact vertices.
%
    
nports = size( port_contacts, 1 );

nedges = size( edges, 1 );    

P = sparse( nedges, nports );

for j=1:nports
    % We will be searching for paths between source and
    % destination vertices to identify the port edges
    allsrc = reshape( port_contacts{ j, 1 }, [], 1 );
    alldst = reshape( port_contacts{ j, 2 }, [], 1 );
    % How many paths?
    npaths = max( length(allsrc), length(alldst) );
    % Start searching!
    src = allsrc;
    dst = alldst;
    for k=1:npaths
        path = shortestpath2( G, 
                              [ src ; allsrc(1:end*isempty(src)) ], ...
                              [ dst ; alldst(1:end*isempty(dst)) ] );
        src( src == path(1)   ) = [];
        dst( dst == path(end) ) = [];

        path_edges        = [ path(1:end-1), path(2:end) ];
        path_edges_sorted = sort( path_edges, 2 );
        edge_sign = -1 + 2*( path_edges(:,1) == path_edges_sorted(:,1) );
        [ found, edge_index ] = ismember(path_edges_sorted, edges, "rows");

        P = P + sparse( edge_index, j, edge_sign ./ npaths, nedges, nports );
    end
end
