function Istiched = myStich(I1, I2)
I1(ismissing(I1)) = 0; % take care of NaN

% padd warped image with zeros to make images in the same size
[h1, w1, ~] = size(I1);
[h2, w2, ~] = size(I2);
I1 = padarray(I1, [h2-h1, w2-w1], 'post');

% find valid pixels in warped image
isValid = I1 > 0;

% % stitching B - if in one channel value is bigger than 0 take all channels,
% % even if their value is 0
% isRedValid = I1(:,:,1) > 0;
% isGreenValid = I1(:,:,2) > 0;
% isBlueValid = I1(:,:,3) > 0;
% temp = logical(isRedValid + isGreenValid + isBlueValid);
% isValid(:,:,1) = temp;
% isValid(:,:,2) = temp;
% isValid(:,:,3) = temp;

% replece pixels with pixels from warped image
Istiched = I2;
Istiched(isValid) = I1(isValid);
