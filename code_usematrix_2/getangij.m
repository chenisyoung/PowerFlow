function angij = getangij(angle_y, alphaU)
% angle_y ���ɾ���ĽǶ�
% aplhaU ��ѹ��� ��Ϊ����
    alphaU_t = alphaU.';
    % ��ʽ i - j - alpha
    angij = repmat(alphaU_t, length(alphaU), 1) - repmat(alphaU, 1, length(alphaU)) - angle_y;
    angij = sparse(angij);
end