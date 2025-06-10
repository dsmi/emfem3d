function [ tri, tri_tetra ] = surftri(p,t)
%SURFTRI Find surface triangles from tetrahedra mesh
%   [ tri, tri_tetra ] = surftri(p,t)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

% Form all faces, non-duplicates are surface triangles
faces=[t(:,[1,2,3]);
       t(:,[1,2,4]);
       t(:,[1,3,4]);
       t(:,[2,3,4])];
node4=[t(:,4);t(:,3);t(:,2);t(:,1)];
tri_tetra=repmat(transpose(1:size(t,1)),4,1);
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
tri=faces(ix(qx),:);
node4=node4(ix(qx));
tri_tetra=tri_tetra(ix(qx));

% Orientation
v1=p(tri(:,2),:)-p(tri(:,1),:);
v2=p(tri(:,3),:)-p(tri(:,1),:);
v3=p(node4,:)-p(tri(:,1),:);
ix=find(dot(cross(v1,v2,2),v3,2)>0);
tri(ix,[2,3])=tri(ix,[3,2]);
