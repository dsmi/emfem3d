function test_tri_u()

%% fd = @(p) -min(min(min(1+p(:,2),1-p(:,2)),1+p(:,1)),1-p(:,1));
%% fh = @(p) ones(size(p,1),1);
%% [ r, tri ] = ...
%% distmesh( fd, fh, 0.5, [-1,-1;1,1], [-1,-1;-1,1;1,-1;1,1] );

% Test triangle
r   = [ 1.2 2.2 ; 1.3 3.5 ; 0.1 2.0 ];
tri = [ 1 2 3 ];

% Test values at the vertices
for t=1:size(tri,1)
    tr = [ r(tri(t,1),:) ; r(tri(t,2),:) ; r(tri(t,3),:) ];
    for j=1:3
        for v=1:3
            u = tri_u( tr, tr(v,:), j );
            ut = ( j == v );
            assert( abs( u - ut ) < 1.0e-10 )
        end
    end
end


%% % Plot last triangle
%% p = simplexquad(20, tr);
%% u = tri_u( tr, p, 3 );

%% colormap("jet");
%% plot( r(1,1), r(1,2), 'k+', r(2,1), r(2,2), 'kx', r(3,1), r(3,2), 'k*' );
%% hold on
%% scatter( p(:,1), p(:,2), [], u );
%% hold off
