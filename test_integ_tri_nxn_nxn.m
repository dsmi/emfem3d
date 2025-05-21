function test_integ_tri_nxn_nxn
    
r = [ 0.2  0.5  0.0 ; ...
      1.7  0    0.3 ; ...
      0.2  1.2  0.4 ; ...
     -0.3  1.4  0.8 ];

tri = [ 1 2 3 ; ...
        2 4 3 ];

for ni = 1:3
    for nj = 1:3
        v  = integ_tri_nxn_nxn( r, tri, ni, nj );
        vt = v*0;
        for t=1:size(tri,1)
            rt    = r( tri(t,:), : );
            vt(t) = iquad_tri_nxn_nxn( rt, ni, nj );
        end
        assert( norm( vt - v ) < norm( vt ) * 1e-10 );
    end
end
