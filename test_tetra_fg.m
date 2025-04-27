function test_tetra_fg()
    
fh = @(p) ones(size(p,1),1);

pfix = [-1,-1,-1;-1,1,-1;1,-1,-1;1,1,-1; -1,-1,1;-1,1,1;1,-1,1;1,1,1];

dppiped = @(p,x1,x2,y1,y2,z1,z2) ...
           -min(min(min(min(min(-y1+p(:,2) ,y2-p(:,2)), ...
                                -x1+p(:,1)),x2-p(:,1)), ...
                                -z1+p(:,3)),z2-p(:,3));

fd = @(p) dppiped(p,-1,1,-1,1,-1,1);

[r,tetra]=distmeshnd(fd,fh,1.0,[-1,-1,-1;1,1,1],pfix);

% tetramesh(tetra,r)

%% r = [ 0, 1, 0 ; ...
%%       0, 0, 0 ; ...
%%       2, 0, 0 ; ...
%%       0, 0, 3 ];

%% tetra = [ 1 2 3 4 ];

vol = tetra_v( r, tetra );

for j=1:1
    for t=1:size(tetra,1)

        % vertices of the current tetrahedron
        rt = r( tetra(t,:), : );
        
        [f, g] = tetra_fg(rt);

        % test points
        p=simplexquad(3,rt);

        [ i1, i2 ] = tetra_edge_verts( j );
        
        [a1,b1,c1,d1] = tetra_abcd( rt, [ 1, 2, 3, 4 ], i1 );
        [a2,b2,c2,d2] = tetra_abcd( rt, [ 1, 2, 3, 4 ], i2 );

        l = tetra_edge_length( rt, [ 1, 2, 3, 4 ], j );

        L1     = 1/(6*vol(t))*(a1 + p(:,1)*b1 + p(:,2)*c1 + p(:,3)*d1);
        gradL1 = 1/(6*vol(t))*[ b1, c1, d1 ]; 

        L2     = 1/(6*vol(t))*(a2 + p(:,1)*b2 + p(:,2)*c2 + p(:,3)*d2);
        gradL2 = 1/(6*vol(t))*[ b2, c2, d2 ];

        np = size(p,1);

        % vector basis functions calculated via abcd
        Na = (repmat(L1,1,3).*repmat(gradL2,np,1) ...
              - repmat(L2,1,3).*repmat(gradL1,np,1)).*repmat(l,np,3);

        % vector basis function calculated via fg
        Nb = repmat(f(j,:),np,1) + cross(repmat(g(j,:),np,1),p,2);

        assert( norm( Na - Nb ) < norm( Na ) * 1e-10 );
    end
end
