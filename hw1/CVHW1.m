% HW1 - Aram Gasparian, 310410865 | Hana Matatov, 203608302
clc, close all;
%% Q2 (there is no Q1 in the PDF)
% section a
axes = -4:0.2:4;
[X,Y] = meshgrid(axes,axes);
Z = X.*exp(-X.^2 - Y.^2);
coordinate = [X(:) Y(:) Z(:)];
figure(1);
surf(X, Y, Z);
title('z = x*exp(-x^2 - y^2)');
xlabel('x');
ylabel('y');
zlabel('z');

% section b
% define 10 coordinates of point of light
S_mat = [4 4 0.2;
    4 -4 0.2;
    -4 4 0.2;
    -4 -4 0.2;
    0 0 0.1;
    0 0 0.2;
    0 2 0.1;
    0 -2 0.1;
    2 0 0.1;
    -2 0 0.1;];

figure(2)
ro = 0.1;
% calculate norms of surface
[Nx, Ny, Nz] = surfnorm(X, Y, Z);
nT = [Nx(:) Ny(:) Nz(:)];

for i = 1:size(S_mat,1)
    S_vec = S_mat(i,:);
    s = (repmat(-S_vec, size(coordinate,1), 1));
    % normalize vector s
    s = s./repmat(sqrt(sum(s.^2,2)), 1,3);
    % calculate irradiance
    Ics = ro*sum(nT.*s, 2);
    I = reshape(Ics, length(axes), []);
    % plot results
    subplot(2,5,i)
    surf(X,Y,I)
    title(['S=(' num2str(S_vec(1)) ',' num2str(S_vec(2)) ',' num2str(S_vec(3)) ')']);
    xlabel('x');
    ylabel('y');
    [U,S,V] = svd(I);
    SVD(i).I = I;
    SVD(i).U = U;
    SVD(i).S = S;
    SVD(i).V = V;
end

% section c
figure(3)
for i = 1:size(S_mat,1)
    subplot(2,5,i)
    eigNum = size(SVD(i).S,2);
    eigValSquare = SVD(i).S(SVD(i).S>0)'.^2;
    rk = [0 cumsum(eigValSquare)/sum(eigValSquare)];
    plot(0:eigNum, rk);
    xlim([0 5]);
    grid on;
    S_vec = S_mat(i,:);
    title(['S=(' num2str(S_vec(1)) ',' num2str(S_vec(2)) ',' num2str(S_vec(3)) ')']);
end

%% Q3 & Q4 are not matlab related

%% Q5
close all; 
% section a
%   calculate image points pts3d & pts2d
inl1;
%   remove NaN
noNaNPoints = ~isnan(pts2d(1,:));
pts2d = pts2d(:,noNaNPoints);
pts3d = pts3d(:,noNaNPoints);
%   build A mat (A*Pcs=0)
A = zeros(length(noNaNPoints));
for j = 1:sum(noNaNPoints)
    X = pts3d(:,j);
    u = pts2d(1,j);
    v = pts2d(2,j);    
    A(2*j-1,:) = [X' zeros(1,4) -X'*u];
    A(2*j,:)   = [zeros(1,4) X' -X'*v];
end

[U,D,V] = svd(A);
P = reshape(V(:,size(D,2)), [4 3])';
Pnorm = P/norm(P); % constraint for optimization problem
% section b
pts2dEst = Pnorm*pts3d;
pts2dEst = pts2dEst./repmat(pts2dEst(3,:), 3, 1); % normalize w curl
figure(1);
plot(pts2dEst(1,:),pts2dEst(2,:),'ro'); % plot projected poiints using matrix P
for i = 1:size(pts2d,2)
    error(i) = norm(pts2d(:,i) - pts2dEst(:,i));
end
E = sum(error)/length(error); % error units are pixels

% section c
[K, R] = rq(Pnorm(1:3,1:3));
scaling = K(3,3);
K = K/scaling; % K(3,3) = 1
R = scaling*R;

% section f
t = K\P(1:3,4);
c = -inv(R)*t;

% section g
figure(2);
plot3(c(1),c(2),c(3),'rx'); % plot camera location in 3D


%% Q6
trainingPath = [pwd '\leaf-data\training\'];
testPath = [pwd '\leaf-data\test\leaf6.png'];
testImg = imread(testPath);
testImgGray = rgb2gray(testImg); % convert target image to gray
sizeVec = size(testImgGray);
% convert training images to gray
for i=1:5
    trainingSet{i} = imread([trainingPath 'leaf' num2str(i) '.png']);
    trainingSetGray{i} = rgb2gray(trainingSet{i});
    sizeVec(i+1,:) = size(trainingSetGray{i});
end
resize = [max(sizeVec(:,1)), max(sizeVec(:,2)) ]; % maximal width and height
Threshold = 200/255;
% imshow(testImgGray);
testBinary = im2bw(testImgGray, Threshold);
testBinary = ~testBinary;
% padd image with zeros
[testSize(1), testSize(2)] = size(testBinary);
testBinary = [zeros(testSize(1),floor((resize(2)-testSize(2))/2)) ,testBinary , zeros(testSize(1),ceil((resize(2)-testSize(2))/2))];
[testSize(1), testSize(2)] = size(testBinary);
testBinary = [zeros(floor((resize(1)-testSize(1))/2),testSize(2)); testBinary; zeros(ceil((resize(1)-testSize(1))/2),testSize(2))];

testBinaryMorph = imclose(testBinary, strel('disk', 30));
scale = sum(sum(testBinaryMorph));
% figure(2)
% imshow(testBinaryMorph);
% figure(3)
% diff = testBinaryMorph-testBinary;
% imshow(diff);
for i=1:5
    trainingSetBinary{i} = ~im2bw(trainingSetGray{i}, Threshold);
    [rowSize, colSize] = size(trainingSetBinary{i});
    trainingSetBinary{i} = [zeros(rowSize,floor((resize(2)-colSize)/2)) ,trainingSetBinary{i} , zeros(rowSize,ceil((resize(2)-colSize)/2))];
    [rowSize, colSize] = size(trainingSetBinary{i});
    trainingSetBinary{i} = [zeros(floor((resize(1)-rowSize)/2),colSize); trainingSetBinary{i}; zeros(ceil((resize(1)-rowSize)/2),colSize)];

    diff = testBinaryMorph - trainingSetBinary{i};
%     figure(i*10);
%     imshow(diff);

% scoring alogirithm, found number of white, grey, black pixels
    blank = (diff == 0); 
    score1 = sum(sum(blank.*testBinaryMorph));
    score0 = sum(sum(diff == 1)); 
    scoreMinus = sum(sum(diff == -1)); 
    score(i) = (score1-scoreMinus-score0)/scale;
    if score(i) < 0
        score(i) = 0;
    end
end
score
[val,idxMatching] = max(score);
fprintf("The matching test leaf is: leaf%d \nWith score of %f\n", idxMatching, val);

