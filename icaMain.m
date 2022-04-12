clear; close all; clc;

dims = 2;
xVecArr = zeros(dims, 10);
mu = 0.001;
el = 30;
phiFn = @(x) tanh(x);

xCount = size(xVecArr, 2);
iMat = eye(dims);
wMatArr = zeros(dims, dims, el);
wMatArr(:, :, 1) = eye(dims);
for i = 1:el - 1
    eMat = zeros(dims, dims);
    for j = 1:xCount
        yVec = wMatArr(:, :, i) * xVecArr(:, j);
        eMat = eMat + phiFn(yVec) * yVec';
    end
    eMat = eMat / xCount;
    wMatArr(:, :, i + 1) = wMatArr(:, :, i) - mu * (eMat - iMat) * wMatArr(:, :, i);
end

yMat = zeros(dims, xCount);
for i = 1:xCount
    yMat(:, i) = wMatArr(:, :, end) * xVecArr(:, i);
end

yTen = zeros(dims, dims, xCount);
for i = 1:xCount
    for j = 1:dims
        yTen(:, j, i) = wMatArr(:, :, end) \ (yMat(:, i) .* iMat(:, j));
    end
end
