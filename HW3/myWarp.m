% section 5 q1 - Image warping
% I - input image
% H - affine transform
function Ip = myWarp(I, H)
[row, col, ~] = size(I);

[x, y] = meshgrid(1:col, 1:row);
origCoordinates = [ x(:) y(:) ones(length(x(:)),1)]';

% find warped image size
prjCorners = H*origCoordinates;
prjCol = ceil(max(prjCorners(1,:)));
prjRow = ceil(max(prjCorners(2,:)));

% create warped image coordinates
[x1, y1] = meshgrid(1:prjCol, 1:prjRow);
prjCoordinates = [ x1(:) y1(:) ones(length(x1(:)),1)]';

% inverse transform
invCoordinates = H\prjCoordinates;
xq = reshape(invCoordinates(1,:), prjRow, prjCol);
yq = reshape(invCoordinates(2,:), prjRow, prjCol);

% interpolate pixel values
Ip = zeros(prjRow, prjCol, 3);
Ip(:,:,1) = interp2(x, y, I(:,:,1), xq, yq, 'bilinear', 0);
Ip(:,:,2) = interp2(x, y, I(:,:,2), xq, yq, 'bilinear', 0);
Ip(:,:,3) = interp2(x, y, I(:,:,3), xq, yq, 'bilinear', 0);
