function v = integ_tri_nxn_nxn( r, tri, ni, nj )
% function v = integ_tri_nxn_nxn( r, tri, ni, nj )
%
%  Evaluates integral of dot product of cross products of the vector basis
%  functions whith the trinagle normal over the triangle:
%   v = \integ ( n \times N_i ) \cdot ( n \times N_i ) dS
%  The normal is calculated as cross( (r1-r3), (r2-r3), 2 ) / length
%  (although because it is  on the both sides of the dot product the sign
%  does not depend on the direction) and the vector basis functions are
%  those defined on a tetrahedron that this triangle is a face of.
%  Notice that the cross product being integrated evaluates to
%    n \times N(p)_i = (v3-p)*l/(2*A)
%  where p is the evaluation point and v3 is the 'free vertex' of the
%  triangle which is one that does not belong to the edge.
%

% Triangle vertices
    r1 = r(tri(:, 1), :);
r2 = r(tri(:, 2), :);
r3 = r(tri(:, 3), :);
    
% Free vertex i and j (one which does not belong to edge i/j)
ri = r(tri(:, ni), :);
rj = r(tri(:, nj), :);

% Triangle areas
nn = cross(r2-r1,r3-r1,2);
A = sqrt(sum(nn.*nn, 2))/2;

% Integral of r over the surface of the triangle, to be multiplied by 2*A
% as 1/6 is the integral of (either of) the barycentric coordinate
integ_r = (r1/6+r2/6+r3/6);

% Integral of dot(r,r) over the triangle, to be multiplied by 2*A
integ_rr = dot(r1,r1,2)/12 + dot(r1,r2,2)/24 + dot(r1,r3,2)/24 + ...
           dot(r2,r1,2)/24 + dot(r2,r2,2)/12 + dot(r2,r3,2)/24 + ...
           dot(r3,r1,2)/24 + dot(r3,r2,2)/24 + dot(r3,r3,2)/12;

li = tri_edge_length( r, tri, ni );
lj = tri_edge_length( r, tri, nj );

% 2*A canceled from numerator and denominator. dot(ri,rj,2) is constant
% and needs to be multiplied by just A to get the integral over the
% triangle, thus it is divided by 2
v = (dot(ri,rj,2)/2-dot(ri+rj,integ_r,2)+integ_rr).*li.*lj./(2*A);
