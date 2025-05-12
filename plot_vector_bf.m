
r = [ 0.1, 1, 0.2 ; ...
      -0.1, 0.3, -0.2 ; ...
      1, 0.5, -0.1 ; ...
      0.5, 0.4, 1 ];

tetra = [ 1 2 3 4 ];

% Edge of the basis function that we are plotting
j=1;

% vertices of our only tetrahedron
rt = r( tetra(1,:), : );
        
[f, g] = tetra_fg(rt);

%% % test points
%% p=simplexquad(3,rt);

qN = 3;
[ qA qW ] = simplexquad(qN, 2);
qA = [ qA 1-sum(qA,2) ]; % barycentric coordinates

v1 = rt(1,:);
v2 = rt(2,:);
v3 = rt(4,:);

p = v1.*qA(:,1) + v2.*qA(:,2) + v3.*qA(:,3);

np = size(p,1);

% vector basis function calculated via fg
Nj = repmat(f(j,:),np,1) + cross(repmat(g(j,:),np,1),p,2);

% Normal and area
n = cross( (v1-v3), (v2-v3), 2 ); % non-normalized
A = sqrt( sum ( n.*n, 2 ) ) / 2 % triangle area
n = repmat( n ./ (2*A), size(Nj, 1), 1 ); % normalized

% Length of the edge
l = sqrt( sum( (v1-v2).*(v1-v2), 2 ) )

nxnj = cross( n, Nj, 2 );


[ i1, i2 ] = tetra_edge_verts( j );

%% tetramesh(tetra,r)

faces=[tetra(:,[1,2,3]);
       tetra(:,[1,2,4]);
       tetra(:,[1,3,4]);
       tetra(:,[2,3,4])];

%% trimesh( faces, r(:,1), r(:,2), r(:,3) )

scatter3( rt([ i1, i2 ],1), rt([ i1, i2 ],2), rt([ i1, i2 ],3), 'b')
hold on
%% scatter3( p(:,1), p(:,2), p(:,3), 'c')
scatter3( [ v1(:,1) v2(:,1) v3(:,1)], ...
          [ v1(:,2) v2(:,2) v3(:,2)], ...
          [ v1(:,3) v2(:,3) v3(:,3)],'c')

quiver3( p(:,1), p(:,2), p(:,3), Nj(:,1), Nj(:,2), Nj(:,3), 'r' );

quiver3( p(:,1), p(:,2), p(:,3), n(:,1), n(:,2), n(:,3), 'm' );

quiver3( p(:,1), p(:,2), p(:,3), nxnj(:,1), nxnj(:,2), nxnj(:,3), 'b' );

% plot edges
for j=1:6
    [ i1, i2 ] = tetra_edge_verts( j );
    plot3( rt([ i1, i2 ],1), rt([ i1, i2 ],2), rt([ i1, i2 ],3), '-g')
end

xlabel('x')
ylabel('y')
zlabel('z')

hold off
