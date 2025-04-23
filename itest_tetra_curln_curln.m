function v = itest_tetra_curln_curln( r, ni, nj )
% function v = itest_tetra_curln_curln( r, ni, nj )
%
%  Evaluates the same integral as integ_tetra_curln_curln:
%   v = \integ \delta \cross N_i \delta \cross N_j dV
%  but using a different approach, to be used for testing.
%  This function can only process one tetrahedron at once,
%   r is 4-by-2 matrix with the vertex coordinates.    
%

% Vertices of our only tetrahedron
tetra = [ 1 2 3 4 ];

vol = tetra_v( r, tetra );

% Edge length
l = zeros( 6, 1 );

% Unit vector of the edge
e = zeros( 6, 3 );

for i=1:6
    l( i ) = tetra_edge_length( r, tetra, i );
    [ i1, i2 ] = tetra_edge_verts( i );
    r1 = r( tetra(:,i1), : );    
    r2 = r( tetra(:,i2), : );
    e(i,:) = ( r2 - r1 ) / l(i);
end

[ f, g ] = tetra_fg( r );

% TFEMEM, alternative expression after (8.68)
v = 4*dot(g(ni,:), g(nj,:), 2)*vol;
