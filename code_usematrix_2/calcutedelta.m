function [unb, maxun] = calcutedelta(n, rundata, nowP, nowQ, Sb, balance, pvs)
% rundata ���в���
% nowP ��ǰ�������������й�
% nowQ ��ǰ�������������޹�
% Sb  ��׼����

P=sparse(rundata(:,1),1,(rundata(:,2)-rundata(:,4))/Sb,n,1);                         
Q=sparse(rundata(:,1),1,(rundata(:,3)-rundata(:,5))/Sb,n,1);%�������Ǽ���PQ�ڵ�ע�빦�ʵģ���Ҫ���Ի�׼����SB

deltaP = P - nowP;
deltaQ = Q - nowQ;
% ��ƽ��ڵ��pv�ڵ���Ӧ������Ϊ0 
deltaP(balance) = 0;
deltaQ(balance) = 0;
deltaQ(pvs) = 0;
unb = [deltaP;deltaQ];
maxun = max(abs(unb));
unb = sparse(unb);
end