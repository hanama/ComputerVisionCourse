%% Q6
trainingPath = 'leaf-data\training\';
testPath = 'leaf-data\test\leaf6.png';
testImg = imread(testPath);
testImgGray = rgb2gray(testImg);
sizeVec = size(testImgGray);
for i=1:5
    trainingSet{i} = imread([trainingPath 'leaf' num2str(i) '.png']);
    trainingSetGray{i} = rgb2gray(trainingSet{i});
    sizeVec(i+1,:) = size(trainingSetGray{i});
end
resize = [max(sizeVec(:,1)), max(sizeVec(:,2)) ];
Threshold = 200/255;
% imshow(testImgGray);
testBinary = im2bw(testImgGray, Threshold);
testBinary = ~testBinary;
[testSize(1), testSize(2)] = size(testBinary);
testBinary = [zeros(testSize(1),floor((resize(2)-testSize(2))/2)) ,testBinary , zeros(testSize(1),ceil((resize(2)-testSize(2))/2))];
[testSize(1), testSize(2)] = size(testBinary);
testBinary = [zeros(floor((resize(1)-testSize(1))/2),testSize(2)); testBinary; zeros(ceil((resize(1)-testSize(1))/2),testSize(2))];
%imshow(I2);

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
    figure(i*10);
    imshow(diff);
    blank = (diff == 0); % need to change
    score1 = sum(sum(blank.*testBinaryMorph));
    score0 = sum(sum(diff == 1)); 
    scoreMinus = sum(sum(diff == -1)); % need to change 
    score(i) = (score1-scoreMinus-score0)/scale;
    if score(i) < 0
        score(i) = 0;
    end
end
score
[val,idxMatching] = max(score);
fprintf("The matching test leaf is: leaf%d \nWith score of %f\n", idxMatching, val);
