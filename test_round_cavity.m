function test_round_cavity

% Run simulation, compare results against golden
    
round_cavity(0)

[ f, S ] = tsread( "round_cavity.z1p" );

[ f0, S0 ] = tsread( "round_cavity_golden.z1p" );

abstol = 1e-2;
reltol = 2e-1;

cmpres = cmpnetwp( f0, S0, f, S, abstol, reltol );
assert( 0 == cmpres );
