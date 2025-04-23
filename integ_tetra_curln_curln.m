function v = integ_tetra_curln_curln( r, tetra, ni, nj )
% function v = integ_tetra_curln_curln( r, tetra, ni, nj )
%
%  Evaluates integral of product of curls of the (vector) basis
%  functions associated with edges i and j in a tetrahedron
%   v = \integ \delta \cross N_i \delta \cross N_j dV
%

vol = tetra_v( r, tetra );

[ i1, i2 ] = tetra_edge_verts( ni );
[ j1, j2 ] = tetra_edge_verts( nj );
    
[ ai1, bi1, ci1, di1 ] = tetra_abcd( r, tetra, i1 );
[ ai2, bi2, ci2, di2 ] = tetra_abcd( r, tetra, i2 );

[ aj1, bj1, cj1, dj1 ] = tetra_abcd( r, tetra, j1 );
[ aj2, bj2, cj2, dj2 ] = tetra_abcd( r, tetra, j2 );

li = tetra_edge_length( r, tetra, ni );
lj = tetra_edge_length( r, tetra, nj );

v = li.*lj./(324*vol.*vol.*vol).* ...
    (   (ci1.*di2 - di1.*ci2).*(cj1.*dj2 - dj1.*cj2) ...
      + (di1.*bi2 - bi1.*di2).*(dj1.*bj2 - bj1.*dj2) ...
      + (bi1.*ci2 - ci1.*bi2).*(bj1.*cj2 - cj1.*bj2) );
