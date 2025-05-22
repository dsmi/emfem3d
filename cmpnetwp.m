function rc = cmpnetwp( freq1, P1, freq2, P2, abstol, reltol )
% rc = cmpnetwp( freq1, P1, freq2, P2, abstol, reltol )
%
% Compare network parameters with tolerance. Returns zero if
% the parameters are the same.
%    

if ( !isequal( size( freq1 ), size( freq2 ) ) )
    rc = 1;
    return;
end

if ( !isequal( size( P1 ), size( P2 ) ) )
    rc = 2;
    return;
end

nf1 = norm( freq1 );
nf2 = norm( freq2 );

if ( abs( nf1 - nf2 ) > abstol )
    rc = 3;
    return;
end

if ( nf1 - nf2 ) > reltol * max( nf1, nf2 ) )
    rc = 4;
    return;
end

np1 = norm( P1 );
np2 = norm( P2 );

if ( abs( np1 - np2 ) > abstol )
    rc = 5;
    return;
end

if ( np1 - np2 ) > reltol * max( np1, np2 ) )
    rc = 6;
    return;
end

rc = 0;
return;
