function angij = getangij(angle_y, alphaU)
% angle_y 导纳矩阵的角度
% aplhaU 电压相角 均为弧度
    alphaU_t = alphaU.';
    % 公式 i - j - alpha
    angij = repmat(alphaU_t, length(alphaU), 1) - repmat(alphaU, 1, length(alphaU)) - angle_y;
    angij = sparse(angij);
end