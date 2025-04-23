function v = integ_tetra_n_n( r, tetra, ni, nj )
% function v = integ_tetra_n_n( r, tetra, ni, nj )
%
%  Evaluates integral of dot the (vector) basis functions
%  associated with edges i and j in a tetrahedron
%   v = \integ N_i \dot N_j dV
%

vol = tetra_v( r, tetra );

% Incorrect orientation will result in wrong sign of the integral.
assert( min(vol) >= 0 );    

li = tetra_edge_length( r, tetra, ni );
lj = tetra_edge_length( r, tetra, nj );

a = b = c = d = cell( 6, 1 );

for j=1:6
    [ a{j}, b{j}, c{j}, d{j} ] = tetra_abcd( r, tetra, j );
end

f = @( i, j ) b{i}.*b{j} + c{i}.*c{j} + d{i}.*d{j};

switch ( [ min( ni, nj ), max( ni, nj ) ] )
    case [ 1, 1 ]
        v = li.*lj./(360*vol).*( f(2,2) - f(1,2) + f(1,1) );
    case [ 1, 2 ]
        v = li.*lj./(720*vol).*( 2*f(2,3) - f(2,1) - f(1,3) + f(1,1) );
    case [ 1, 3 ]
        v = li.*lj./(720*vol).*( 2*f(2,4) - f(2,1) - f(1,4) + f(1,1) );
    case [ 1, 4 ]
        v = li.*lj./(720*vol).*( f(2,3) - f(2,2) - 2*f(1,3) + f(1,2) );
    case [ 1, 5 ]
        v = li.*lj./(720*vol).*( f(2,2) - f(2,4) - f(1,2) + 2*f(1,4) );
    case [ 1, 6 ]
        v = li.*lj./(720*vol).*( f(2,4) - f(2,3) - f(1,4) + f(1,3) );
    case [ 2, 2 ]
        v = li.*lj./(360*vol).*( f(3,3) - f(1,3) + f(1,1) );
    case [ 2, 3 ]
        v = li.*lj./(720*vol).*( 2*f(3,4) - f(1,3) - f(1,4) + f(1,1) );
    case [ 2, 4 ]
        v = li.*lj./(720*vol).*( f(3,3) - f(2,3) - f(1,3) + 2*f(1,2) );
    case [ 2, 5 ]
        v = li.*lj./(720*vol).*( f(2,3) - f(3,4) - f(1,2) + f(1,4) );
    case [ 2, 6 ]
        v = li.*lj./(720*vol).*( f(1,3) - f(3,3) - 2*f(1,4) + f(3,4) );
    case [ 3, 3 ]
        v = li.*lj./(360*vol).*( f(4,4) - f(1,4) + f(1,1) );
    case [ 3, 4 ]
        v = li.*lj./(720*vol).*( f(3,4) - f(2,4) - f(1,3) + f(1,2) );
    case [ 3, 5 ]
        v = li.*lj./(720*vol).*( f(2,4) - f(4,4) - 2*f(1,2) + f(1,4) );
    case [ 3, 6 ]
        v = li.*lj./(720*vol).*( f(4,4) - f(3,4) - f(1,4) + 2*f(1,3) );
    case [ 4, 4 ]
        v = li.*lj./(360*vol).*( f(3,3) - f(2,3) + f(2,2) );
    case [ 4, 5 ]
        v = li.*lj./(720*vol).*( f(2,3) - 2*f(3,4) - f(2,2) + f(2,4) );
    case [ 4, 6 ]
        v = li.*lj./(720*vol).*( f(3,4) - f(3,3) - 2*f(2,4) + f(2,3) );
    case [ 5, 5 ]
        v = li.*lj./(360*vol).*( f(2,2) - f(2,4) + f(4,4) );
    case [ 5, 6 ]
        v = li.*lj./(720*vol).*( f(2,4) - 2*f(2,3) - f(4,4) + f(3,4) );
    case [ 6, 6 ]
        v = li.*lj./(360*vol).*( f(4,4) - f(3,4) + f(3,3) );
end

