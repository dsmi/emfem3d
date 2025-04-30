
r = [ 0, 1, 0 ; ...
      0, 0, 0 ; ...
      1, 0.5, 0 ; ...
      0.5, 0.5, 1 ];

tetra = [ 1 2 3 4 ];

% Edge of the basis function that we are plotting
j=1;

% vertices of our only tetrahedron
rt = r( tetra(1,:), : );
        
[f, g] = tetra_fg(rt);

% test points
p=simplexquad(3,rt);

np = size(p,1);

% vector basis function calculated via fg
Nj = repmat(f(j,:),np,1) + cross(repmat(g(j,:),np,1),p,2);


[ i1, i2 ] = tetra_edge_verts( j );

%% tetramesh(tetra,r)

faces=[tetra(:,[1,2,3]);
       tetra(:,[1,2,4]);
       tetra(:,[1,3,4]);
       tetra(:,[2,3,4])];

%% trimesh( faces, r(:,1), r(:,2), r(:,3) )

scatter3( rt([ i1, i2 ],1), rt([ i1, i2 ],2), rt([ i1, i2 ],3), 'b')
hold on
%% scatter3( p(:,1), p(:,2), p(:,3), 'r')

quiver3( p(:,1), p(:,2), p(:,3), Nj(:,1), Nj(:,2), Nj(:,3), 'r' );

% plot edges
for j=1:6
    [ i1, i2 ] = tetra_edge_verts( j );
    plot3( rt([ i1, i2 ],1), rt([ i1, i2 ],2), rt([ i1, i2 ],3), '-g')
end

xlabel('x')
ylabel('y')
zlabel('z')

hold off
