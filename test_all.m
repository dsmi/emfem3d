
testfuncs = { ...
             @test_det4, ...
             @test_integ_tetra_curln_curln, ...
             @test_integ_tetra_n_n, ...
             @test_tetra_abcd, ...
             @test_tetra_edge_length, ...
             @test_tetra_fg, ...
             @test_tetra_v, ...
             @test_tri_u, ...
             @test_extrude_mesh, ...
             @test_ptonrect, ...
             @test_collect_tetra_edges, ...
             @test_shortestpath2, ...
             @test_tsread, ...
             @test_tri_edge_length, ...
             @test_integ_tri_nxn_nxn, ...
             @test_collect_tri_edges, ...
             @test_surftri ...
};

startTime = tic;

for tf=testfuncs
    disp( [ 'Running ', func2str( tf{1} ) ] );
    tf{1}();
end

elapsedTime = toc( startTime );

disp( [ 'All done in ', num2str( elapsedTime ), ' sec.' ] );

