clear; close all; clc;

% mu : ステップサイズμ
% el : 繰り返し回数L
mu = 0.5;
el = 30;

% xVecArr : 入力信号xの列
% fs : サンプリング周波数
% xLen : 入力信号長T
% dim : 次元N
% iMat : N次単位行列I
[xVecArr, fs] = audioread("1+2+3.wav");
[xLen, dim] = size(xVecArr);
iMat = eye(dim);

% pFn : 生成モデルp(y)
% phiFn : スコア関数φ(y)
% jFn : カルバック・ライブラー・ダイバージェンスの計算関数J
pFn = @(y) 1 / pi * sech(y);
phiFn = @(y) tanh(y);
jFn = @(w, y) -log(abs(det(w))) - sum(log(pFn(y)), "all") / xLen;

% wMat : 分離行列W
% yVecArr : 分離信号yの列
% jArr : カルバックーライブラ・ダイバージェンスJの列
wMat = iMat;
yVecArr = xVecArr;
jArr = zeros(el, 1);

jArr(1) = jFn(wMat, yVecArr);
for i = 1:el - 1
    eMat = phiFn(yVecArr') * yVecArr / xLen;
    wMat = wMat - mu * (eMat - iMat) * wMat;
    yVecArr = (wMat * xVecArr')';
    jArr(i + 1) = jFn(wMat, yVecArr);
end

plot(jArr);

maxVol = max(yVecArr, [], "all");
yVecArr = yVecArr / maxVol * 0.8;
audiowrite("1&2&3.wav", yVecArr, fs);