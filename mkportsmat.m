function P = mkportsmat( G, edges, port_contacts )

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
