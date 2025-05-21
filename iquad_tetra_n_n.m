function v = iquad_tetra_n_n( r, ni, nj )
% v = iquad_tetra_n_n( r, ni, nj )
%
%  Evaluates the same integral as integ_tetra_n_n
%   v = \integ N_i \dot N_j dV
%  using a quadrature; to be used for testing.
%  This function can only process one tetrahedron at once,
%   r is 4-by-3 matrix with the vertex coordinates.    
%

% quadrature points
[ qp, qw ] = simplexquad( 3, r );

% coefficients to calculate the basis functions
[ f, g ] = tetra_fg( r );

nqp = size(qp,1);

Ni = repmat(f(ni,:),nqp,1) + cross(repmat(g(ni,:),nqp,1),qp,2);

Nj = repmat(f(nj,:),nqp,1) + cross(repmat(g(nj,:),nqp,1),qp,2);

v = sum( dot( Ni, Nj, 2 ).*qw, 1 );
