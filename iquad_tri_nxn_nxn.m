function v = iquad_tri_nxn_nxn( r, ni, nj )
% v = iquad_tri_nxn_nxn( r, ni, nj )
%
%  Evaluates the same integral as integ_tri_nxn_nxn
%   v = \integ ( n \times N_i ) \cdot ( n \times N_i ) dS
%  using a quadrature; to be used for testing.
%  This function can only process one triangle at once,
%   r is 3-by-3 matrix with the vertex coordinates, column
%   is xyz row is vertex.
%

r1 = r(1,:);
r2 = r(2,:);
r3 = r(3,:);

% Quadrature points
qN = 5;
[ qA qW ] = simplexquad( qN, 2 );
qA = [ qA 1-sum(qA,2) ]; % barycentric coordinates

p = r1.*qA(:,1) + r2.*qA(:,2) + r3.*qA(:,3);

np = size(p,1);
    
% Normal and area
nn = cross( (r1-r3), (r2-r3), 2 ); % non-normalized
A = sqrt( sum ( nn.*nn, 2 ) ) / 2; % triangle area
n = nn ./ (2*A);                   % normalized

%% % temp to debug integ_tri_nxn_nxn
%% iquad_r  = sum(p.*qW,1)
%% iquad_rr = sum(dot(p,p,2).*qW,1)

%% % temp -- calculate integral as integ_tri_nxn_nxn to debug it
%% ll = sqrt(sum((r2 - r3).^2,2))
%% rr = r1
%% iquad0 = sum(sum(( rr - p ).*( rr - p ), 2).*qW, 1 )*2*A.*ll./(2*A).*ll./(2*A)
%% iquad1 = sum(sum(( rr.*rr - 2*rr.*p + p.*p ), 2).*qW, 1 )*2*A.*ll./(2*A).*ll./(2*A)
%% k = 2*A.*ll./(2*A).*ll./(2*A);
%% iquad2 = ( sum( rr.*rr, 2 )/2 - dot( 2*rr, iquad_r, 2 ) + iquad_rr )*k

% Make a 'fake' tetrahedron to calculate the basis function
% values over the face
rt = [ r1 - nn; r1; r2; r3 ];

[ f, g ] = tetra_fg( rt );

% calculate vector basis functions
tri2tetra = [ 6, 5, 4 ]; % triangle edge -> edge of our fake tetrahedron
i = tri2tetra( ni );
j = tri2tetra( nj );
Ni = repmat(f(i,:),np,1) + cross(repmat(g(i,:),np,1),p,2);
Nj = repmat(f(j,:),np,1) + cross(repmat(g(j,:),np,1),p,2);

% cross with the normal
nxNi = cross( repmat( n, np, 1 ), Ni, 2 );
nxNj = cross( repmat( n, np, 1 ), Nj, 2 );

% dot and sum the quadrature points
v = sum( dot( nxNi, nxNj, 2 ).*qW*A*2, 1 );
