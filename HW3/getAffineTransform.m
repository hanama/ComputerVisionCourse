% section 3 q1
function H = getAffineTransform(p1, p2)
% P*h=PCurl
% Create P matrix
pp = [p1' ones(size(p1',1),1)];
pX = [pp zeros(size(p1',1),3)];
pX = reshape(pX',3,[])';
pY = [zeros(size(p1',1),3) pp];
pY = reshape(pY',3,[])';
P = [pX pY];
PCurl = p2(:);
% solve using Least Squares method
h = (P'*P)\P'*PCurl;
H = [h(1:3)'; h(4:6)'; 0 0 1];


