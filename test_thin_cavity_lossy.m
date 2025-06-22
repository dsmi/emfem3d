function test_thin_cavity_lossy

% Run simulation, compare results against golden
    
thin_cavity_lossy(0)

[ f, S ] = tsread( "thin_cavity_highloss.z2p" );

[ f0, S0 ] = tsread( "thin_cavity_highloss_golden.z2p" );

abstol = 1e-2;
reltol = 2e-1;

cmpres = cmpnetwp( f0, S0, f, S, abstol, reltol );
assert( 0 == cmpres );
