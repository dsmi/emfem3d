function [ r, tetra ] = extrude_mesh( rt, tri, zl )
% [ r, tetra ] = extrude_mesh( rt, tri, zl )
%
% 'Extrude' 2d mesh in z direction to make 3d one.
% zl is the z coordinates of the extruded layers boundaries.
%

% number of layers minus one (number of the layer iterfaces)
nl = length(zl);

% Number of vertices per layer
nvl = size( rt, 1 );

% tetrahedra mesh vertices
r = [ repmat(rt, nl, 1), ...
      reshape(repmat(reshape(zl, 1, nl), nvl, 1), nl*nvl, 1) ];

% Number of triangles per layer
ntri = size( tri, 1 );

tetra = ones( ntri*(nl-1)*3, 4 );

for l=1:(nl-1)
    
    v1 = tri(:,1) + nvl*(l-1);
    v2 = tri(:,2) + nvl*(l-1);
    v3 = tri(:,3) + nvl*(l-1);
    v4 = tri(:,1) + nvl*l;
    v5 = tri(:,2) + nvl*l;
    v6 = tri(:,3) + nvl*l;

    tetra1 = [ v1, v2, v3, v4 ];
    tetra2 = [ v2, v3, v4, v5 ];
    tetra3 = [ v3, v4, v5, v6 ];
    
    tetra( (ntri*3*(l-1)+0*ntri+1):(ntri*3*(l-1)+1*ntri), : ) = tetra1;
    tetra( (ntri*3*(l-1)+1*ntri+1):(ntri*3*(l-1)+2*ntri), : ) = tetra2;
    tetra( (ntri*3*(l-1)+2*ntri+1):(ntri*3*(l-1)+3*ntri), : ) = tetra3;
end

