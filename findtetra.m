function tidx = findtetra( r, tetra, fd )
% tidx = findtetra( r, tetra, fd )
%
%  Find tetra with centers inside the volume described by
%  the signed distance function fd (ones for which fd is
%  negative)

% Centers
rc = (1/4)*(r(tetra(:,1),:)+r(tetra(:,2),:)+r(tetra(:,3),:)+r(tetra(:,4),:));
    
tidx = find( fd(rc) < 0 );
