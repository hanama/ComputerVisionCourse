%%%%%%%%%% Q2 - Photometry %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%A%%%%%%%%%%%%%%%%%%%%%
[x,y] = meshgrid(-5:0.2:5,-5:0.2:5);
z = x .* exp(-(x.*x)-(y.*y));
surf(x,y,z);

%%%%%%%%%%%%%%%%%%B%%%%%%%%%%%%%%%%%%%%%
[nx, ny, nz] = surfnorm(x,y,z);
n =  reshape([nx ny nz], 51, 51, 3);

r = 0.1;

for k = [1:10]
    
    s = randn(3,1);
    sn = s./norm(s);
    
    res = zeros(51, 51);

    for i = [1:51]
        for j = [1:51]

            np = n(i, j, :);
            np = reshape(np, 1, 3);

            irr = r * np * sn;
            res(i,j) = irr;
        end
    end
    subplot(2,5,k);
    surf(x, y, res);
    %%%%%%%%%%%%%%%%%%C%%%%%%%%%%%%%%%%%%%%%
    [U,S,V] = svd(res);
    rs = rank(S);
    title(arrayfun(@(x) {num2str(x)},sn));
    rs
end 

%%
%%%%%%%%%%Q6%%%%%%%%%%%%
l1 = imread('leaf1.png');
l2 = imread('leaf2.png');
l3 = imread('leaf3.png');
l4 = imread('leaf4.png');
l5 = imread('leaf5.png');
l6 = imread('leaf6.png');

imshow(l6)

lg1 = rgb2gray(l1);
lg2 = rgb2gray(l2);
lg3 = rgb2gray(l3);
lg4 = rgb2gray(l4);
lg5 = rgb2gray(l5);
lg6 = rgb2gray(l6);

lgbw6 = ~im2bw(lg6);

figure(2)
imshow(lgbw6);

se = strel('disk',30);

closeBW = imclose(lgbw6,se);

figure(3)
imshow(closeBW);


lgbw1 = ~im2bw(lg1);
lgbw2 = ~im2bw(lg2);
lgbw3 = ~im2bw(lg3);
lgbw4 = ~im2bw(lg4);
lgbw5 = ~im2bw(lg5);