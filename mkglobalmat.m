function [ K, b, S ] = mkglobalmat( r, tetra, edges, tetrae, tetraes, ...
                                    tri, tri_tetra, trie, tries, ...
                                    w, er, lt, fref, sigma, ...
                                    P, fpec, fsurfz, fmatch )
% [ K, b, S ] = mkglobalmat( r, tetra, edges, tetrae, tetraes, ...
%                            tri, tri_tetra, trie, tries, ...
%                            w, er, lt, fref, sigma, ...
%                            P, fpec, fsurfz, fmatch )
%
% Makes the finite elements global/system matrix.
%  r, tetra, edges, tetrae, tetraes -- tetrahedral mesh
%  tri, tri_tetra, trie, tries      -- surface triangles
%  w, er, lt, fref                  -- parameters of dielectric, per tetrahedron
%  sigma                            -- conductivity for surface impedance b.c.
%  P                                -- ports matrix
% fpec, fsurfz, fmatch              -- distance functions for pec, surface impedance
%                                      and maching boundary b.c. correspondingly

nedges = size( edges, 1 );

% length of the edges
edgelen = sqrt( sum( (r(edges(:,2),:) - r(edges(:,1),:)).^2, 2) );

% Unknown (non-pec) edges. Check if both ends and center of the
% edge is on the pec boundary to see if this is a pec edge. 
eunk = find( ~( abs( fpec( (r(edges(:,1),:)+r(edges(:,2),:))*.5) ) < 1e-10 ...
                & abs( fpec( r(edges(:,2),:) ) ) < 1e-10 ...
                & abs( fpec( r(edges(:,1),:) ) ) < 1e-10 ) );

% number of the 'unknown' edges
neunk = length(eunk);

% Matrix to get rows/columns for the unknown edges from full K
S = sparse( eunk, transpose( 1:neunk ), ones( neunk, 1 ), nedges, neunk );

% properties of the medium
epsd = eps0*debye(er, lt, 2*pi*fref, w);
Z0 = sqrt(mu0./epsd);
k0 = w*sqrt(epsd*mu0);

% the same over the entire domain
k0_Z0 = w*mu0;

% Full admittance matrix K
K = sparse( nedges, nedges );

% Tetrahedra
for ni=1:6
    for nj=1:6
        K0 = integ_tetra_curln_curln( r, tetra, ni, nj );
        K2 = integ_tetra_n_n( r, tetra, ni, nj );
        Ke = tetraes(:,ni).*tetraes(:,nj).*( K0 - k0.*k0.*K2 );
        K = K + sparse( tetrae(:,ni), tetrae(:,nj), Ke, nedges, nedges );
    end
end

% Identify triangles for surface impedance or matching boundary
% condition to be applied.

vsurfz = find( abs( fsurfz( r ) ) < 1e-7 );

tsurfz = find( ismember( tri(:,1), vsurfz ) ...
             & ismember( tri(:,2), vsurfz ) ...
             & ismember( tri(:,3), vsurfz ) );

vmatch = find( abs( fmatch( r ) ) < 1e-7 );

tmatch = find( ismember( tri(:,1), vmatch ) ...
             & ismember( tri(:,2), vmatch ) ...
             & ismember( tri(:,3), vmatch ) );

% Surface impedance or matching boundary or nothing as defined
% by fsurfz and fmatch
ktri = zeros(size(tri,1),1);
ktri(tsurfz) = sqrt(j*sigma*k0_Z0);      % surface impedance
ktri(tmatch) = j*k0(tri_tetra(tmatch));  % matching
for ni=1:3
    for nj=1:3
        Ks = ktri.*integ_tri_nxn_nxn( r, tri, ni, nj );
        Ke = tries(:,ni).*tries(:,nj).*Ks;
        K = K + sparse( trie(:,ni), trie(:,nj), Ke, nedges, nedges );
    end
end

% This drops known pec edges from the global matrix
K = S'*K*S;

% Rhs vector
b = full( S'*spdiags( -j*k0_Z0*edgelen, 0, nedges, nedges )*P );
