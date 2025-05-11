function drawports(P,edges,r)

nports = size(P,2);
                   
for k=1:nports
    [pe,pj,pw] = find(P(:,k));
    r0 = r(edges(pe,1),:);
    rr = r(edges(pe,2),:) - r(edges(pe,1),:);
    nwidx = find( pw < 0 ); % ones with neg weight are to be flipped
    r0(nwidx,:) = r(edges(pe(nwidx),2),:);
    rr(nwidx,:) = r(edges(pe(nwidx),1),:) - r(edges(pe(nwidx),2),:);
    quiver3( r0(:,1), r0(:,2), r0(:,3), rr(:,1), rr(:,2), rr(:,3), 0 );
end
