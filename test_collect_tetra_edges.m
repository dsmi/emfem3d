function test_collect_tetra_edges

% Test mesh
fh = @(p) ones(size(p,1),1);

pfix = [-1,-1,-1;-1,1,-1;1,-1,-1;1,1,-1; -1,-1,1;-1,1,1;1,-1,1;1,1,1];

dppiped = @(p,x1,x2,y1,y2,z1,z2) ...
           -min(min(min(min(min(-y1+p(:,2) ,y2-p(:,2)), ...
                                -x1+p(:,1)),x2-p(:,1)), ...
                                -z1+p(:,3)),z2-p(:,3));

fd = @(p) dblock(p,-1,1,-1,1,-1,1);

[r,tetra]=distmeshnd(fd,fh,0.5,[-1,-1,-1;1,1,1],pfix);

%% tetramesh(tetra,r)

[ edges, tetrae, tetraes ] = collect_tetra_edges( tetra );

assert( sort(edges,2) == edges )
assert( unique(edges,'rows') == edges )

% Validate trie and tries
for t=1:size(tetra,1)
    for j=1:6
        [ i1, i2 ] = tetra_edge_verts( j );
        v1 = tetra(t,i1);
        v2 = tetra(t,i2);
        assert( 1 == tetraes(t,j) || -1 == tetraes(t,j) );
        if 1 == tetraes(t,j)
            assert( [ v1, v2 ] == edges( tetrae(t,j ), : ) );
        else
            assert( [ v2, v1 ] == edges( tetrae(t,j ), : ) );
        end
    end
end
