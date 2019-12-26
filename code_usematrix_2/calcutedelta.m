function [unb, maxun] = calcutedelta(rundata, nowP, nowQ, Sb, balance, pvs)
% rundata ���в���
% nowP ��ǰ�������������й�
% nowQ ��ǰ�������������޹�
% Sb  ��׼����

P = rundata(:,2) - rundata(:,4); % Pi
P = P ./ Sb;
Q = rundata(:,3) - rundata(:,5); % Qi
Q = Q ./ Sb;

deltaP = P - nowP;
deltaQ = Q - nowQ;
% ��ƽ��ڵ��pv�ڵ���Ӧ������Ϊ0
deltaP(balance) = 0;
deltaQ(balance) = 0;
deltaQ(pvs) = 0;
unb = [deltaP;deltaQ];
maxun = max(abs(unb));

%     unbalance_P=P-real(S);           %ȡ������ʵ�������й���ƽ����
%     unbalance_Q=Q-imag(S);           %ȡ�������鲿�����޹���ƽ����
%     unbalance_P(balance)=0;           %��ƽ��ڵ���й���ƽ������Ϊ��
%     unbalance_P=sparse(unbalance_P);
%     unbalance_Q(balance)=0;           %��ƽ��ڵ���޹���ƽ������Ϊ��
%     unbalance_Q=sparse(unbalance_Q);
%     unbalance_Q(pvs)=0;           %��PV�ڵ���޹���ƽ������Ϊ��
%     unbalance_Q=sparse(unbalance_Q);
%     unbalance=sparse([unbalance_P;unbalance_Q]);         %���й��޹���ƽ������Ϊһ������,Ϊ������ȡ��������׼��
%     maxunb=max(abs(unbalance));%�ж����ƽ��������Ϊֹͣ����������
end