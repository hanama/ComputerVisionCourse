% HW3 - Task2
%Hana Matatov, 203608302, shanama@campus.technion.ac.il
%Aram Gasparian, 310410865, aram89@campus.technion.ac.i
clc, close all;
%------------------------------------------------------------
addpath('sift', 'Stop Images');
%% Q1 - SIFT
% read images 
files = dir('Stop Images\*.jpg');
imageList = {files(:).name};
I = cell(4);
descr = cell(4);
frames = cell(4);
% section 1 - extract sift descriptor
for i = 1:numel(imageList)
   imageName = imageList{i};
   I{i} = im2double(imread(imageName));
   [frames{i}, descr{i}] = sift(rgb2gray(I{i}));
end

% section 2 - find matching key-points
match12 = siftmatch(descr{1}, descr{2});
match13 = siftmatch(descr{1}, descr{3});
match14 = siftmatch(descr{1}, descr{4});
%%
figure, plotmatches(I{1},I{2}, frames{1}, frames{2}, match12);
figure, plotmatches(I{1},I{3}, frames{1}, frames{3}, match13);
figure, plotmatches(I{1},I{4}, frames{1}, frames{4}, match14);
%%
[H12, inliers12, ransacMatch12] = getTransform(frames{1}, frames{2}, match12);
figure, plotmatches(I{1},I{2}, frames{1}, frames{2}, ransacMatch12);

[H13, inliers13, ransacMatch13] = getTransform(frames{1}, frames{3}, match13);
figure, plotmatches(I{1},I{3}, frames{1}, frames{3}, ransacMatch13);

[H14, inliers14, ransacMatch14] = getTransform(frames{1}, frames{4}, match14);
figure, plotmatches(I{1},I{4}, frames{1}, frames{4}, ransacMatch14);

disp(['Inliers12 = ' num2str(inliers12)]);
disp(H12);
disp('---------------------------------------');
disp(['Inliers13 = ' num2str(inliers13)]);
disp(H13);
disp('---------------------------------------');
disp(['Inliers14 = ' num2str(inliers14)]);
disp(H14);
disp('---------------------------------------');

%%
warpI12 = myWarp(I{1}, H12);
figure, imshow(warpI12);
stiched12 = myStich(warpI12, I{2});
figure, imshow(stiched12);

warpI13 = myWarp(I{1}, H13);
figure, imshow(warpI13);
stiched13 = myStich(warpI13, I{3});
figure, imshow(stiched13);


warpI14 = myWarp(I{1}, H14);
figure, imshow(warpI14);
stiched14 = myStich(warpI14, I{4});
figure, imshow(stiched14);


%% Q2 - Optical flow
N = 16;
T = 1;
% frameIdx = 1;
seq = im2double(imread('seq.gif'));
[h, w, ch, t] = size(seq);
% show all frames as one image
seq = squeeze(seq);
seq2 = reshape(seq,h,[]);
% calculate deriratives for x,y,t
seqDx = seq2 - [seq2(:,2:end) seq2(:,1)];
seqDy = seq2 - [seq2(2:end,:); seq2(1,:)];
seqDt = seq2 - [seq2(:,w+1:end) seq2(:,1:w)];
% Ix/Iy size:[patches in frame, size of patch, frames]
patches = im2col(seqDx, [N N], 'distinct');
Ix = reshape(patches, N*N, h*w/N^2, t);

patches = im2col(seqDy, [N N], 'distinct');
Iy = reshape(patches, N*N, h*w/N^2, t);

patches = im2col(seqDt, [N N], 'distinct');
It = reshape(patches, N*N, h*w/N^2, t);

frameNum = 10;
Frames = [];
iptsetpref('ImshowBorder','tight');
for frameIdx = 1:frameNum
    % iterate over regions in one frame
    [u, patchIdx] = calcU(h, w, N, T, Ix, Iy, It, frameIdx);
    h1 = figure('visible','off'); imshow(seq(:,:,frameIdx));
    hold all;
    % calculate top left pixel coordinates of patches
    y = (floor((patchIdx-1)./(h/N))').*N+1;
    x = mod((patchIdx-1)', h/N).*N+1;
    quiver(y,x,u(1,:)',u(2,:)', 'Color', 'red');
%     title(['N = ' num2str(N) ', T = ' num2str(T)]);
    hold off;
    currFrame = getframe(h1);
    Frames = [Frames currFrame.cdata];
end
figure, imshow(Frames);
%%
% section 3 - T = 0.1
T = 0.1;
% iterate over regions in one frame
[u, patchIdx] = calcU(h, w, N, T, Ix, Iy, It, frameIdx);
figure, imshow(seq(:,:,frameIdx));
hold all;
% calculate top left pixel coordinates of patches 
y = (floor((patchIdx-1)./(h/N))').*N+1;
x = mod((patchIdx-1)', h/N).*N+1;
quiver(y,x,u(1,:)',u(2,:)', 'Color', 'red');
title(['N = ' num2str(N) ', T = ' num2str(T)]);
hold off;

% section 4 - N = 8, T = [1, 0.1]
N = 8;
patches = im2col(seqDx, [N N], 'distinct');
Ix = reshape(patches, N*N, h*w/N^2, t);
 
patches = im2col(seqDy, [N N], 'distinct');
Iy = reshape(patches, N*N, h*w/N^2, t);

patches = im2col(seqDt, [N N], 'distinct');
It = reshape(patches, N*N, h*w/N^2, t);

Tvec = [0.1 1];
for T = Tvec
    % iterate over regions in one frame
    [u, patchIdx] = calcU(h, w, N, T, Ix, Iy, It, frameIdx);    
    figure, imshow(seq(:,:,frameIdx));
    hold all;
    % calculate top left pixel coordinates of patches
    y = (floor((patchIdx-1)./(h/N))').*N+1;
    x = mod((patchIdx-1)', h/N).*N+1;
    quiver(y,x,u(1,:)',u(2,:)', 'Color', 'red');
    title(['N = ' num2str(N) ', T = ' num2str(T)]);
end

%%