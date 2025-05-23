%% function test_collect_tetra_edges

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

edges = collect_tetra_edges( tetra );

tri = surftri(r,tetra);

[ trie, tries ] = collect_tri_edges( tri, edges );

% Validate trie and tries
for t=1:size(tri,1)
    for j=1:3
        [ i1, i2 ] = tri_edge_verts( j );
        v1 = tri(t,i1);
        v2 = tri(t,i2);
        assert( 1 == tries(t,j) || -1 == tries(t,j) );
        if 1 == tries(t,j)
            assert( [ v1, v2 ] == edges( trie(t,j ), : ) );
        else
            assert( [ v2, v1 ] == edges( trie(t,j ), : ) );
        end
    end
end
