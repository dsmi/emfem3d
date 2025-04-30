function [ sp ] = shortestpath2( G, s, d )
% [ sp ] = shortestpath2( G, s, d )
%
% Finds shortest path between (multiple) s and d in a weighted graph
% described by adjacency matrix G, which can be sparse. G(i,j) is length
% of the edge from j to i (opposite to the usual convention to make the
% matrix access more efficient), zero value means no path between nodes.
% Edge lengths must be positive.
% Uses a simplistic implementation of the Dijkstra's algorithm without
% a priority queue which works in O(v^2) time (thus 2 in the name, which
% however should be ok if the sought paths are not too long.
% sp is the path vertices.

N = size(G,1);
    
done = zeros(N, 1);      % 'queue', 1 if we finised examining the vertex
dist = inf*ones(N, 1);   % shortes distance to any of the sources
prev = zeros(N, 1);      % predecessor, 0 for all the sources

% 1 if the vertex is a destination
is_dest = zeros(N, 1 );
is_dest(d) = 1;

% Start the search!
dist(s) = 0;

while sum(done)~=N
    unseen_dist = dist;
    unseen_dist(find(done)) = inf;
    [ udist u ] = min(unseen_dist);
    done(u) = 1;
    if is_dest(u)
        reachd = u;
        done(:) = 1;
    else
        [ni,nj,nw] = find(G(:,u));
        for k=1:length(ni)
            v = ni(k);
            altdist = dist(u)+nw(k);
            if altdist < dist(v)
                dist(v)=altdist;
                prev(v)=u;
            end
        end
    end
end

% Reverse iteration to build the path
sp = [ reachd ];
while prev(sp(1)) ~= 0
    sp = [ prev(sp(1)) ; sp ];
end
