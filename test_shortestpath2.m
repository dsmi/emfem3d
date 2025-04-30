function test_shortestpath2

% Test graph:
%      (1)     (1)     (1)
%  v1------v2------v3------v4
%          |       |
%          |(2)    |(10)
%     (7)  |   (3) |   (5)
%  v5------v6------v7------v8
%          |       |        |
%          |(2)    |(3)     |(3)
%          |   (1) |    (5) |
%          v9------v10-----v11

e = [ 1 2; 2 3; 3 4; 2 6; 3 7; 5 6; 6 7; 7 8; 6 9; 7 10; 8 11; 9 10; 10 11 ];
w = [   1;   1;   1;   2;  10;   7;   3;   5;   2;    3;    3;    1;     5 ];

nedges = size( e, 1 );

% Adjacency matrix as shortestpath2 wants it
G = sparse( [ e(:,1), e(:,2) ], ...
            [ e(:,2), e(:,1) ], ...
            [      w,    w   ], nedges, nedges );

assert( [ 5; 6 ]              == shortestpath2( G, 5, 6 ) );
assert( [ 8; 7; 6; 2; 3; 4  ] == shortestpath2( G, 8, 4 ) );
assert( [ 2; 6; 7; 8  ]       == shortestpath2( G, [ 2 3 ], [ 8 11 ] ) );
