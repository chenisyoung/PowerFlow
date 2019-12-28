clc;clear;
[sysdata, linedata, branchdata, transferdata, rundata, pvdata]=loadfile();
    % һЩ���õĲ���
    e = sysdata(2, 1);
    maxdiedai = sysdata(1, 4);
    nodes = sysdata(1,1); % �ڵ���
    Sb = sysdata(1,3); %���ʱ���ֵ ��׼
    balance = sysdata(3,2); % ƽ��ڵ���
    pvs = pvdata(:,1); % pv�ڵ���
    Y=Y_matrix(sysdata, linedata, branchdata, transferdata); % ʹ�ø���ֱ�Ӽ��� �ٶȲ�� �����׿���
    % Y = Form_Y_matrix(nodes, linedata(:,2:end), transferdata(:,2:end), branchdata);%�������ʹ��ָ�����ϵ��㷨

    %% ���ֵ ֮��Ū��function
    tic
    U = ones(sysdata(1,1),1);       % ����ֵ1
    alphaU = zeros(sysdata(1,1),1);  % ����λ0
    U(pvdata(:,1))=pvdata(:,2); % pv�ڵ��ѹ������ 
    U_polar = U.*exp(1i.*alphaU); 
    % abs_y = abs(Y);
    %% ��ʼ�������
    loops = 0;
    while 1
        % �ȼ��� ��ǰ����
        [nowP, nowQ] = calculatepower(U_polar, Y); 
        % ��ƽ���� �� ��ǰ����
        [unb, maxun] = calcutedelta(nodes, rundata, nowP, nowQ, Sb, balance, pvs);
        
        if maxun < e  % ��������㹻
            disp('���ȴﵽҪ��');
            disp(['ѭ������' ': ' num2str(loops)]);
            alphaU = rad2deg(alphaU);
            break;
        end
        if loops > maxdiedai
            disp('�޷�����')
            disp(['ѭ������: ' num2str(loops)]);
            break;
        end
        Jac = Jacobin(nodes, U_polar, Y, balance, pvs);  % �ο�����д��
        % Jac = jacobi(abs_y, angij, U, nodes, balance, pvs); % ����ʽ��д Ч���е��
        % 
        delta = Jac \ unb;
        alphaU = alphaU - delta(1:nodes);  % �������
        U = U - delta(nodes+1:end) .* abs(U);  % ��ֵ
        U_polar = sparse(U.*exp(1i.*alphaU));  % ���»�Ϊ����
        loops = loops + 1;
    end
    toc