function [ f, g ] = tetra_fg(r)
% [ f, g ] = tetra_fg(r)
%
%  Calculate coefficients which implement yet another way to calculate
%  a linear basis function in a tetrahedron:
%   N_j = f_j + cross( g_j, r )
%  Where f_j and g_j are vectors calculated by this function, and r is
%  the point where we want to evaluate the basis function (assumed to be
%  within the tetrahedron, otherwise the basis function is zero).
%  Note that the second coefficient, g, is g(j,:) = curl(N_j)/2 
%  
%  Used for testing.
%
    
% Vertices of our only tetrahedron
tetra = [ 1 2 3 4 ];

vol = tetra_v( r, tetra );

% Edge length
l = zeros( 6, 1 );

% Unit vector of the edge
e = zeros( 6, 3 );

% Edge endpoints
r1 = zeros( 6, 3 );
r2 = zeros( 6, 3 );

for i=1:6
    l( i ) = tetra_edge_length( r, tetra, i );
    [ i1, i2 ] = tetra_edge_verts( i );
    r1(i,:) = r( tetra(:,i1), : );    
    r2(i,:) = r( tetra(:,i2), : );
    e(i,:) = ( r2(i,:) - r1(i,:) ) ./ l( i );
end

f = zeros( 6, 3 );
g = zeros( 6, 3 );

% formulas (8.60) from TFEMEM:
for i=1:6    
    f(7-i,:) = l(7-i)/(6*vol)*cross( r1(i,:), r2(i,:), 2 );
    g(7-i,:) = l(i)*l(7-i)*e(i,:)/(6*vol);
end
