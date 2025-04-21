function [ u gradu ] = tri_u( r, p, j )
% function [ u gradu ] = tri_u( r, p, j )
%
%  Calculates scalar basis function for the specified triangle.
%   r - 3-by-2 -- vertices of the triangle
%   p - N-by-2 -- points where to evaluate the function
%   j - basis function number which is the vertex that the function is 1
%
%  The basis functions returned by this functions are also the area or
%  barycentric coordinates in the triangle where i'th basis function is
%  the weight of the corresponding triangle vertex.
%

cycind = @( i ) rem( i - 1, 3 ) + 1;

x = r(:,1);
y = r(:,2);

a = b = c = zeros( 3, 1 );
for i = 1:3
    a(i) = x(cycind(i+1))*y(cycind(i+2)) - y(cycind(i+1))*x(cycind(i+2));
    b(i) = y(cycind(i+1)) - y(cycind(i+2));
    c(i) = x(cycind(i+2)) - x(cycind(i+1));
end

A = 1/2*(b(1)*c(2) - b(2)*c(1));

u = 1/(2*A)*(a(j) + b(j)*p(:,1) + c(j)*p(:,2));

gradu = 1/(2*A)*[ b(j) c(j) ];
