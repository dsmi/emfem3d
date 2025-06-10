function test_surftri

% Test mesh
fh = @(p) ones(size(p,1),1);

pfix = [-1,-1,-1;-1,1,-1;1,-1,-1;1,1,-1; -1,-1,1;-1,1,1;1,-1,1;1,1,1];

fd = @(p) dblock(p,-1,1,-1,1,-1,1);

[r,tetra]=distmeshnd(fd,fh,0.5,[-1,-1,-1;1,1,1],pfix);

[ tri, tri_tetra ] = surftri( r, tetra );

% make sure that triangle concides with one of the faces of the
% its source tetra by comparing the centers.
for t=1:size(tri,1)

    tt = tri_tetra(t);
    faces=[tetra(tt,[1,2,3]);
           tetra(tt,[1,2,4]);
           tetra(tt,[1,3,4]);
           tetra(tt,[2,3,4])];
    
    c = (r(tri(t,1),:)+r(tri(t,2),:)+r(tri(t,3),:))*(1/3);
    
    fc = (r(faces(:,1),:)+r(faces(:,2),:)+r(faces(:,3),:))*(1/3);

    d = sum(abs(fc - repmat(c,4,1)),2);

    assert( d(1) < 1e-10 || d(2) < 1e-10 || d(3) < 1e-10 || d(4) < 1e-10 );

end
