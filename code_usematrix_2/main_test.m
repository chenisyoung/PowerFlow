clc
clear
[sysdata, linedata, branchdata, transferdata, rundata, pvdata, gendata]=loadfile();
    Y=Y_matrix(sysdata, linedata, branchdata, transferdata);
    % һЩ���õĲ���
    e = sysdata(2, 1);
    maxdiedai = sysdata(1, 4);
    nodes = sysdata(1,1); % �ڵ���
    Sb = sysdata(1,3); %���ʱ���ֵ ��׼
    balance = sysdata(3,2); % ƽ��ڵ���
    pvs = pvdata(:,1); % pv�ڵ���
    abs_y = abs(Y); % ���ɷ�ֵ
    angle_y = angle(Y); % �������
    %% ���ֵ
    U = ones(sysdata(1,1),1);       % ����ֵ1
    alphaU = zeros(sysdata(1,1),1);  % ����λ0
    U(pvdata(:,1))=pvdata(:,2); % pv�ڵ��ѹ������ 
    %% ��ʼ�������
    loops = 1;
    while 1
        angij = getangij(angle_y, alphaU);  % ��ȡ i - j - y
        % ��Ӧ�ȼ��� ��ǰ����
        % [nowP, nowQ] = calculatepower(abs_y, nodes, alphaU, U, balance, pvs); 
        [nowP, nowQ] = calculatepower(alphaU, U, Y, balance, pvs); 
        % ��ƽ���� �� ��ǰ����
        [unb, maxun] = calcutedelta(rundata, nowP, nowQ, Sb, balance, pvs);
        if maxun < e  % ��������㹻
            disp(['���ȴﵽҪ��' ' ' num2str(loops)]);
            break;
        end
        if loops > maxdiedai
            disp(['�޷�����' ' ' num2str(loops)]);
            break;
        end
        Jac = jacobi(abs_y, angij, U, nodes, balance, pvs);  % ���ſɱȾ���

        % 
        delta = Jac \ unb;
        alphaU = alphaU - delta(1:nodes);  % �������
        U = U - delta(nodes+1:end);  % ��ֵ
        loops = loops + 1;
    end