function [ a, b, c, d ] = tetra_abcd(r,tetra,j)
% [ a, b, c, d ] = tetra_abcd(r,tetra,j)
%
%  Calculate coefficients for a linear basis function on a tetrahedron
%  which can be used to calculate the basis function value as
%   u = a + x*b + y*c + z*d
%
%  The idea is to solve the 4-by-4 system enforcing the basis function values
%  at the vertices using Cramer's rule
%
    
x1 = r(tetra(:,1),1);
x2 = r(tetra(:,2),1);
x3 = r(tetra(:,3),1);
x4 = r(tetra(:,4),1);

y1 = r(tetra(:,1),2);
y2 = r(tetra(:,2),2);
y3 = r(tetra(:,3),2);
y4 = r(tetra(:,4),2);

z1 = r(tetra(:,1),3);
z2 = r(tetra(:,2),3);
z3 = r(tetra(:,3),3);
z4 = r(tetra(:,4),3);

I = x1*0+1;

% Volume of the tetrahedron
V = det4(  I,  I,  I,  I, ...
          x1, x2, x3, x4, ...
          y1, y2, y3, y4, ...
          z1, z2, z3, z4 ) / 6;

% Values of the basis functions at the corresponding vertex
u1 = x1*0 + ( 1 == j );
u2 = x1*0 + ( 2 == j );
u3 = x1*0 + ( 3 == j );
u4 = x1*0 + ( 4 == j );

a = det4( u1, u2, u3, u4, ...
          x1, x2, x3, x4, ...
          y1, y2, y3, y4, ...
          z1, z2, z3, z4 ) ./ (6 * V);

b = det4(  I,  I,  I,  I, ...
          u1, u2, u3, u4, ...
          y1, y2, y3, y4, ...
          z1, z2, z3, z4 ) ./ (6 * V);

c = det4(  I,  I,  I,  I, ...
          x1, x2, x3, x4, ...
          u1, u2, u3, u4, ...
          z1, z2, z3, z4 ) ./ (6 * V);

d = det4(  I,  I,  I,  I, ...
          x1, x2, x3, x4, ...
          y1, y2, y3, y4, ...
          u1, u2, u3, u4 ) ./ (6 * V);
