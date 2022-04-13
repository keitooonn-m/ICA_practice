clear; close all; clc;

% dims : 信号源・観測点数N
% xVecArr : 入力信号x
% mu : ステップサイズμ
% el : 繰り返し回数L
% pFn : 生成モデルp
% phiFn : スコア関数φ
dims = 2;
xVecArr = zeros(dims, 10);
mu = 0.001;
el = 30;
pFn = @(y) 1 / pi : sech(y);
phiFn = @(y) tanh(y);

% xLen : 入力信号長
% iMat : N次単位行列
% yVecArr : 分離信号y
% wMatArr : 分離行列の列
xLen = size(xVecArr, 2);
iMat = eye(dims);
yVecArr = zeros(dims, xLen);
wMatArr = zeros(dims, dims, el);

wMatArr(:, :, 1) = iMat;
for i = 1:el - 1
    % eMat : 期待値の行列E
    eMat = zeros(dims, dims);

    % y = W * x
    yVecArr = wMatArr(:, :, i) * xVecArr;
    % E = SUM(R) / T
    for j = 1:xLen
        yVec = wMatArr(:, :, i) * xVecArr(:, j);
        rMat = phiFn(yVec) * yVec';
        eMat = eMat + rMat;
    end
    eMat = eMat / xLen;

    % W[i + 1] = W[i] - μ * (E - I) * W[i]
    wMatArr(:, :, i + 1) = wMatArr(:, :, i) - mu * (eMat - iMat) * wMatArr(:, :, i);
end

% カルバックーライブラ・ダイバージェンスJの列
jArr = zeros(xLen, 1);
for i = 1:el
    jArr(i) = -log(abs(det(wMatArr(:, :, i)))) - sum(log(pFn(yVecArr))) / xLen;
end

plot(jArr);

for i = 1:xLen
    yVecArr(:, i) = wMatArr(:, :, end) * xVecArr(:, i);
end
