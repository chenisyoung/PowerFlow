clc
clear
[sysdata, linedata, branchdata, transferdata, rundata, pvdata, gendata]=loadfile();
tic
    Y=Y_matrix(sysdata, linedata, branchdata, transferdata);
    % ���ʱ���ֵ ��׼
    Sb = sysdata(1,3);
    fuzhi = abs(Y);
    jiao = angle(Y);
    %% ���ֵ
    U = ones(sysdata(1,1),1);       % ����ֵ1
    alphaU = zeros(sysdata(1,1),1);  % ����λ0
    U(pvdata(:,1))=pvdata(:,2); % pv�ڵ��ѹ������ 
    %% ��ʼ�������
    %for i = 1:3 % ��������
        [ddelta, dU] = Jacobian(fuzhi, jiao, sysdata, rundata, U, alphaU, pvdata, Sb);
        U = U - dU.*U;
        alphaU = alphaU - ddelta;
        % disp(['�ڼ��ε���' ' ' num2str(i)]);
    %end
toc