function test_coupled_microstrip

% Run the coupled microstrip simulation, compare results against golden
    
coupled_microstrip(0);

[ f, S ] = tsread( "coupled_microstrip.s4p" );

[ f0, S0 ] = tsread( "coupled_microstrip_golden.s4p" );

abstol = 1e-2;
reltol = 1e-2;

assert( 0 == cmpnetwp( f0, S0, f, S, abstol, reltol ) );

