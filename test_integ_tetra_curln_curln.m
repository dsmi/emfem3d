function test_integ_tetra_curln_curln()
    
fh = @(p) ones(size(p,1),1);

pfix = [-1,-1,-1;-1,1,-1;1,-1,-1;1,1,-1; -1,-1,1;-1,1,1;1,-1,1;1,1,1];

dppiped = @(p,x1,x2,y1,y2,z1,z2) ...
           -min(min(min(min(min(-y1+p(:,2) ,y2-p(:,2)), ...
                                -x1+p(:,1)),x2-p(:,1)), ...
                                -z1+p(:,3)),z2-p(:,3));

fd = @(p) dppiped(p,-1,1,-1,1,-1,1);

[r,tetra]=distmeshnd(fd,fh,0.5,[-1,-1,-1;1,1,1],pfix);

% tetramesh(tetra,r)

%% tetra = tetra(1,:);

%% r = [ 0, 1, 0 ; ...
%%       0, 0, 0 ; ...
%%       2, 0, 0 ; ...
%%       0, 0, 3 ];

%% tetra = [ 1 2 3 4 ];

for ni=1:6
    for nj=1:6
        v = integ_tetra_curln_curln( r, tetra, ni, nj );
        vt = v*0;
        for t=1:size(tetra,1)
            rt = r( tetra(t,:), : );
            vt(t) = itest_tetra_curln_curln( rt, ni, nj );
        end
        assert( norm( vt - v ) < norm( vt ) * 1e-10 );
    end
end
