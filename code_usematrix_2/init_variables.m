%%__��ƽ��������������ʼֵ__
%%__����:����__
%%__�������:2019��12��15��__
%%����˵��
%% n:�ڵ���      pvdata:pv�ڵ�����
function [U_polar, U, alphaU] = init_variables(n, pvdata)
    U = ones(n,1);       % ����ֵ1
    alphaU = zeros(n,1);  % ����λ0
    U(pvdata(:,1))=pvdata(:,2); % pv�ڵ��ѹ������ 
    U_polar = U.*exp(1i.*alphaU); 
end