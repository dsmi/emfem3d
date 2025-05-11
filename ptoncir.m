function p = ptoncir( cx, cy, r, n )
% p = ptoncir( cx, cy, r, n )
%
% Generates n points which lie on the circumference of a circle.
%

a = linspace( 0, 2*pi*(n-1)/n, n-1 );

p = [ cx + r*cos(a'), cy + r*sin(a') ];
