clear; close all; clc;

% mu : ステップサイズμ
% el : 繰り返し回数L
% pFn : 生成モデルp(y)
% phiFn : スコア関数φ(y)
mu = 0.8;
el = 30;
pFn = @(y) 1 / pi * sech(y);
phiFn = @(y) tanh(y);

% xVecArr : 入力信号xの列
% xLen : 入力信号長T
% dim : 次元N
% iMat : N次単位行列I
[xVecArr, fs] = audioread("1+2+3.wav");
xLen = size(xVecArr, 1);
dim = size(xVecArr, 2);
iMat = eye(dim);

% 分離行列Wの列
wMatArr = zeros(dim, dim, el);
wMatArr(:, :, 1) = iMat;
for i = 1:el - 1
    eMat = zeros(dim);
    for j = 1:xLen
        yVec = wMatArr(:, :, i) * xVecArr(j, :)';
        pVec = phiFn(yVec);
        rMat = pVec * yVec';
        eMat = eMat + rMat;
    end
    eMat = eMat / xLen;

    % W[i + 1] = W[i] - μ * (E - I) * W[i]
    wMatArr(:, :, i + 1) = wMatArr(:, :, i) - mu * (eMat - iMat) * wMatArr(:, :, i);
end

% y = W * x
yVecArr = zeros(xLen, dim);
for i = 1:xLen
    yVecArr(i, :) = wMatArr(:, :, end) * xVecArr(i, :)';
end

% カルバックーライブラ・ダイバージェンスJの列
jArr = zeros(el, 1);
for i = 1:el
    jArr(i) = -log(abs(det(wMatArr(:, :, i)))) - sum(log(pFn(yVecArr)), "all") / xLen;
end
plot(jArr);

maxVol = max(yVecArr, [], "all");
yVecArr = yVecArr / maxVol * 0.8;
audiowrite("1&2&3.wav", yVecArr, fs);