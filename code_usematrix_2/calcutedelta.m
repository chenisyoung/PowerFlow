%%__���㲻ƽ����__
%%__����:����__
%%__�������:2019��12��15��__
%%����˵��
%% n:�ڵ���        rundata:���в���,�����ڵ㹦��     
%% nowP:��ǰ��������ĸ��ڵ��й�����
%% nowQ:��ǰ��������ĸ��ڵ��޹�����
%% Sb:��׼����      balance:ƽ��ڵ���      pvs:ƽ��ڵ���
%% unb:���ڵ㲻ƽ����      maxun:���ƽ����
function [unb, maxun] = calcutedelta(n, rundata, nowP, nowQ, Sb, balance, pvs)
% �����и�Ҫע��ĵط� 14�ڵ��30�ڵ��PQ���ǰ��ڵ���Ÿ��� 400�ĺ�1047�Ĵ�����
% ֮ǰû��ȡ��ų��� ֱ��Ĭ�ϰ����,�ʳ��� ���˺ܾ�bug���ҵ�
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