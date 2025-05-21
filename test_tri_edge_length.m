function test_tri_edge_length

r = [ 0, 3, 0 ; ...
      4, 0, 0 ; ...
      0, 0, 0 ];

tri = [ 1 2 3 ];

l1 = tri_edge_length( r, tri, 1 );
assert( abs( l1 - 4 ) < 1.0e-10 );

l2 = tri_edge_length( r, tri, 2 );
assert( abs( l2 - 3 ) < 1.0e-10 );

l3 = tri_edge_length( r, tri, 3 );
assert( abs( l3 - 5 ) < 1.0e-10 );
