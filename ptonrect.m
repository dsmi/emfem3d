function p = ptonrect( cx, cy, dx, dy, nx, ny, spf )
% p = ptonrect( cx, cy, dx, dy, nx, ny, spf )
%
% Generates points which lie on the perimeter of a rectangle.
% Dimensions of the rectangle are (-dx/2 - dx/2), (-dy/2 - dy/2)
% along x and y correspondingly; center ix cx, cy.
% Total number of points generated is 2*nx + 2*ny, with (nx+1) and
% (ny+1) being the number of points on the x and y directed sides.
% If spf is not specified the points on each edge are spaced with
% linspace, spf is used otherwise. It is expected to accept the same
% arguments as linspace and return the row vector of values.
%

if ~exist('spf')
    spf = @linspace;
end

p1 = [ cx + spf(-dx/2, dx/2, nx+1)' cy + repmat(-dy/2, nx+1, 1) ];

p2 = [ cx + repmat(dx/2, ny+1, 1)   cy + spf(-dy/2, dy/2, ny+1)' ];

p3 = [ cx + spf(dx/2, -dx/2, nx+1)' cy + repmat(dy/2, nx+1, 1) ];

p4 = [ cx + repmat(-dx/2, ny+1, 1)  cy + spf(dy/2, -dy/2, ny+1)' ];

p = [ p1(1:end-1,:) ; p2(1:end-1,:); p3(1:end-1,:); p4(1:end-1,:) ];
