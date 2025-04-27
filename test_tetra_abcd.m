function test_tetra_abcd()
    
fh = @(p) ones(size(p,1),1);

pfix = [-1,-1,-1;-1,1,-1;1,-1,-1;1,1,-1; -1,-1,1;-1,1,1;1,-1,1;1,1,1];

dppiped = @(p,x1,x2,y1,y2,z1,z2) ...
           -min(min(min(min(min(-y1+p(:,2) ,y2-p(:,2)), ...
                                -x1+p(:,1)),x2-p(:,1)), ...
                                -z1+p(:,3)),z2-p(:,3));

fd = @(p) dppiped(p,-1,1,-1,1,-1,1);

[r,tetra]=distmeshnd(fd,fh,1.0,[-1,-1,-1;1,1,1],pfix);

% tetramesh(tetra,r)

vol = tetra_v(r, tetra);

for j=1:4
    [a,b,c,d] = tetra_abcd(r,tetra,j);
    for t=1:size(tetra,1)
        x=r(tetra(t,:),1);
        y=r(tetra(t,:),2);
        z=r(tetra(t,:),3);
        A = [ ones( 4, 1 ), x, y, z ];
        abcd_all = A\(6*vol(t)*eye(4,4));
        abcd = abcd_all(:,j);
        assert( abs( a(t) - abcd(1) ) < abs( abcd(1) ) * 1e-10 + 1.0e-10 );
        assert( abs( b(t) - abcd(2) ) < abs( abcd(2) ) * 1e-10 + 1.0e-10 );
        assert( abs( c(t) - abcd(3) ) < abs( abcd(3) ) * 1e-10 + 1.0e-10 );
        assert( abs( d(t) - abcd(4) ) < abs( abcd(3) ) * 1e-10 + 1.0e-10 );
    end
end
