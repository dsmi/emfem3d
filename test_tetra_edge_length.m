%% function test_tetra_edge_length

%%    v1 *
%%     / | \_
%%    |  |   \_
%%    |  |     \_
%%    |  |       \_
%%   |   |         \_
%%   |   |           \_     Y ^ 
%%   |   |             \      |
%%  |    *--------------*     +--> X
%%  |     v2    ___---  v3   /
%%  |   ___ ----            Z
%%  *---                  
%%  v4

r = [ 0, 1, 0 ; ...
      0, 0, 0 ; ...
      2, 0, 0 ; ...
      0, 0, 3 ];

tetra = [ 1 2 3 4 ];

% v1 - v2
l1 = tetra_edge_length( r, tetra, 1 );
assert( abs( l1 - 1 ) < 1.0e-10 );

% v1 - v3
l2 = tetra_edge_length( r, tetra, 2 );
assert( abs( l2 - sqrt( 1^2 + 2^2 ) ) < 1.0e-10 );

% v1 - v4
l3 = tetra_edge_length( r, tetra, 3 );
assert( abs( l3 - sqrt( 1^2 + 3^2 ) ) < 1.0e-10 );

% v2 - v3
l4 = tetra_edge_length( r, tetra, 4 );
assert( abs( l4 - sqrt( 0^2 + 2^2 ) ) < 1.0e-10 );

% v4 - v2
l5 = tetra_edge_length( r, tetra, 5 );
assert( abs( l5 - sqrt( 0^2 + 3^2 ) ) < 1.0e-10 );

% v3 - v4
l6 = tetra_edge_length( r, tetra, 6 );
assert( abs( l6 - sqrt( 2^2 + 3^2 ) ) < 1.0e-10 );
