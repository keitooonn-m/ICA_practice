clear; close all; clc;

dimn = 2;
xMat = zeros(dimn, 10);
mu = 0.001;
el = 30;
phiFn = @(x) tanh(x);

xCount = size(xMat, 2);
iMat = eye(dimn);
wTen = zeros(dimn, dimn, el);
wTen(:, :, 1) = eye(dimn);
for i = 1:el - 1
    eMat = zeros(dimn, dimn);
    for j = 1:xCount
        yVec = wTen(:, :, i) * xMat(:, j);
        eMat = eMat + phiFn(yVec) * yVec';
    end
    eMat = eMat / xCount;
    wTen(:, :, i + 1) = wTen(:, :, i) - mu * (eMat - iMat) * wTen(:, :, i);
end

yMat = zeros(dimn, xCount);
for i = 1:xCount
    yMat(:, i) = wTen(:, :, end) * xMat(:, i);
end

yTen = zeros(dimn, dimn, xCount);
for i = 1:xCount
    for j = 1:dimn
        yTen(:, j, i) = wTen(:, :, end) \ (yMat(:, i) .* iMat(:, j));
    end
end
