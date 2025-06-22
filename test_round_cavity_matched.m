function test_round_cavity_matched

% Run simulation, compare results against golden
    
round_cavity_matched(0)

[ f, S ] = tsread( "round_cavity_matched.z1p" );

[ f0, S0 ] = tsread( "round_cavity_matched_golden.z1p" );

abstol = 1e-2;
reltol = 2e-1;

cmpres = cmpnetwp( f0, S0, f, S, abstol, reltol );
assert( 0 == cmpres );
