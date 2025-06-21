function s = z2s(z)
% s = z2s(z)
%
% Convert Z network parameters to S
% s and z are n-by-n-by-nfreqs
%
% S = inv(Z+I) * (Z-I)

I = diag(ones(1, size(z,1)));

for i=1:size(z,3)
    s(:,:,i) = inv(z(:,:,i)+I) * (z(:,:,i)-I);
end
