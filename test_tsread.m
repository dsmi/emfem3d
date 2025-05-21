function test_tsread

[ fz, Z, tz ] = tsread( 'test_touchstone.z2p' );

assert( strcmp( 'z', tz ) );
assert( 133 == length( fz ) );
assert( [ 2 2 133 ] == size( Z ) );

assert(abs(1.0e+9-fz(1))<1e-5)
assert(abs(1.072267222010323e+09-fz(3))<1e-5)
assert(abs(1.0e+11-fz(133))<1e-1)

assert(abs((2.979445997483481e-07-j*2.320071797647215e+03) - Z(1,1,1))<1e-10)
assert(abs((-1.360611959461646e-07-j*2.242542833518569e+03) - Z(1,2,2))<1e-10)
assert(abs((-1.662668455604388e-06+j*1.121733581060866e+02) - Z(2,1,133))<1e-10)
assert(abs((5.133559919787952e-06-j*5.580472603966135e+01) - Z(2,2,133))<1e-10)

[ fs, S, ts ] = tsread( 'test_touchstone.s4p' );

assert( strcmp( 's', ts ) );
assert( 9 == length( fs ) );
assert( [ 4 4 9 ] == size( S ) );

s12 = S(1,2,3);
assert( abs( 6.299587822380916e-08+j*6.378678054213011e-08 - s12 ) < 1e-12 );





