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

if ( abs( nf1 - nf2 ) > abstol + reltol * max( nf1, nf2 ) )
    rc = 3;
    return;
end

nfreq = size(P1,3);

for i=1:nfreq

    np1 = norm( P1(:,:,i) );
    np2 = norm( P2(:,:,i) );

    if ( abs( np1 - np2 ) > abstol + reltol * max( np1, np2 ) )
        rc = 4;
        return;
    end
    
end

rc = 0;
return;
