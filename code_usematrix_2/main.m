%%__ţ������ѷ�㷨�ĳ�������__
%%__����:����__
%%__�������:2020��1��2��__
%% ����˵��
%%n:�ڵ���        linedata:��·����       branchdata:֧·����
%%transferdata:��ѹ������
clc;clear;
%% ��ȡ����
[sysdata, linedata, branchdata, transferdata, rundata, pvdata]=loadfile();
    % ȡ��һЩ���õĲ���
    e = sysdata(2, 1); % ��������
    MaxIteration = sysdata(1, 4); % ����������
    nodes = sysdata(1,1); % �ڵ���
    Sb = sysdata(1,3); %���ʱ���ֵ ��׼
    balance = sysdata(3,2); % ƽ��ڵ���
    pvs = pvdata(:,1); % pv�ڵ���
tic % ��ʼ��ʱ
%% �γ��ſɱȾ��󲢸���������ֵ
    Y=Y_matrix(nodes, linedata, branchdata, transferdata); % ʹ�ø���ֱ�Ӽ���
    [U_polar, U, alphaU] = init_variables(nodes, pvdata);  % ��ʼ��ֵ
%% ��ʼ�������
    loops = 0; % ��¼ѭ������
    while 1
        [nowP, nowQ] = calculatepower(U_polar, Y); % �ȼ��� ��ǰ���ڵ㹦�� �� Pi Qi
        % ��ƽ���� �� ��ǰ����(���ƽ����)
        [unb, maxun] = calcutedelta(nodes, rundata, nowP, nowQ, Sb, balance, pvs);
        % loss(loops+1) = maxun; % �������ƽ�����½�����
        if maxun < e  % ��������㹻
            disp('���ȴﵽҪ��');
            disp('���ȴﵽҪ��ʱ��');
            toc % ��¼����ʱ��
            disp(['ѭ������' ': ' num2str(loops)]);
            outputresult(nowP, nowQ, U, alphaU, U_polar, Y, linedata);
            alphaU = rad2deg(alphaU); % ������ת��Ϊ�Ƕ�
            break;
        end
        if loops > MaxIteration  % ��������������ʱ������޷�����
            disp('�˴μ����޷������޷�����')
            disp(['ѭ������: ' num2str(loops)]);
            break;
        end
        Jac = Jacobin(nodes, U_polar, Y, balance, pvs);  % �����ſɱȾ���
        delta = Jac \ unb;      % �� inv(Jac) * unb 
        alphaU = alphaU - delta(1:nodes);  % ��ǵ�����
        U = U - delta(nodes+1:end) .* abs(U);  % ��ֵ������
        U_polar = sparse(U.*exp(1i.*alphaU));  % ���»�Ϊ������ϡ��
        loops = loops + 1;
    end
    disp('����ͻ�ͼ��ɺ�ʱ��');
    toc % ��¼�������Լ���ͼ֮���ʱ��