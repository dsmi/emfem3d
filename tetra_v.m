function v = tetra_v(r,tetra)
% v = tetra_abcd(r,tetra)
%
%  Calculate volume of the tetrahedrons.
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
v = det4(  I,  I,  I,  I, ...
          x1, x2, x3, x4, ...
          y1, y2, y3, y4, ...
          z1, z2, z3, z4 ) / 6;
