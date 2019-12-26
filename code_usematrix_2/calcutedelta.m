function [unb, maxun] = calcutedelta(rundata, nowP, nowQ, Sb, balance, pvs)
% rundata 运行参数
% nowP 当前参数计算所得有功
% nowQ 当前参数计算所得无功
% Sb  基准容量

P = rundata(:,2) - rundata(:,4); % Pi
P = P ./ Sb;
Q = rundata(:,3) - rundata(:,5); % Qi
Q = Q ./ Sb;

deltaP = P - nowP;
deltaQ = Q - nowQ;
% 将平衡节点和pv节点相应的量设为0
deltaP(balance) = 0;
deltaQ(balance) = 0;
deltaQ(pvs) = 0;
unb = [deltaP;deltaQ];
maxun = max(abs(unb));

%     unbalance_P=P-real(S);           %取复功率实部计算有功不平衡量
%     unbalance_Q=Q-imag(S);           %取复功率虚部计算无功不平衡量
%     unbalance_P(balance)=0;           %将平衡节点的有功不平衡量置为零
%     unbalance_P=sparse(unbalance_P);
%     unbalance_Q(balance)=0;           %将平衡节点的无功不平衡量置为零
%     unbalance_Q=sparse(unbalance_Q);
%     unbalance_Q(pvs)=0;           %将PV节点的无功不平衡量置为零
%     unbalance_Q=sparse(unbalance_Q);
%     unbalance=sparse([unbalance_P;unbalance_Q]);         %将有功无功不平衡量合为一个矩阵,为后面求取修正量做准备
%     maxunb=max(abs(unbalance));%判断最大不平衡量，作为停止迭代的条件
end