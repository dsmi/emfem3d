function test_dipole_antenna

% Run simulation, compare results against golden
    
dipole_antenna(0)

[ f, S ] = tsread( "dipole_antenna.z1p" );

[ f0, S0 ] = tsread( "dipole_antenna_golden.z1p" );

abstol = 1e-2;
reltol = 1e-1;

cmpres = cmpnetwp( f0, S0, f, S, abstol, reltol );
assert( 0 == cmpres );
