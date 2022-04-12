function outTen = ica(inMat, mu, el, phiFn)
    % inMat : 入力行列
    % mu : ステップサイズμ
    % el : 繰り返し回数L
    % phiFn : スコア関数

    % dims : 信号源・観測点数N
    % inLen : 信号長
    % iMat : N次単位行列
    dims = size(inMat, 1);
    inLen = size(inMat, 2);
    iMat = eye(dims);

    % wTen : 分離行列の列
    wTen = zeros(dims, dims, el);
    wTen(:, :, 1) = eye(dims);
    for i = 1:el - 1
        % eMat : 標本期待値の行列
        eMat = zeros(dims, dims);
        for j = 1:inLen
            % NxN * Nx1 => Nx1
            yVec = wTen(:, :, i) * inMat(:, j);
            % NxN + Nx1 * 1xN => NxN
            eMat = eMat + phiFn(yVec) * yVec';
        end
        eMat = eMat / inLen;
        
        wTen(:, :, i + 1) = wTen(:, :, i) - mu * (eMat - iMat) * wTen(:, :, i);
    end
    
    yMat = zeros(dims, inLen);
    for i = 1:inLen
        yMat(:, i) = wTen(:, :, end) * inMat(:, i);
    end

    outTen = zeros(dims, dims, inLen);
    for i = 1:inLen
        for j = 1:dims
            outTen(:, j, i) = wTen(:, :, end) \ (yMat(:, i) .* iMat(:, j));
        end
    end
end
