% For question 2
function [u, patchIdx] = calcU(h, w, N, T, Ix, Iy, It, frameIdx)
u = [];
patchIdx = [];
% iterate over regions in one frame
for i=1:h*w/N^2
    % compute A matrix and b vector   (Au=b)
   A = [Ix(:,i,frameIdx) Iy(:,i,frameIdx)];
   b = -It(:,i,frameIdx);
   eigenValues = eig(A'*A);
   if min(eigenValues) > T % check if the minimum eigen value is bigger then threshold
       u(:,end+1) = (A'*A)\(A'*b);
       patchIdx(end+1) = i; % save patch indexes of calculated u vector
   end
end