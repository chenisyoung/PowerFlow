function [unb, maxun] = calcutedelta(n, rundata, nowP, nowQ, Sb, balance, pvs)
% rundata 运行参数
% nowP 当前参数计算所得有功
% nowQ 当前参数计算所得无功
% Sb  基准容量

P=sparse(rundata(:,1),1,(rundata(:,2)-rundata(:,4))/Sb,n,1);                         
Q=sparse(rundata(:,1),1,(rundata(:,3)-rundata(:,5))/Sb,n,1);%这两步是计算PQ节点注入功率的，都要除以基准容量SB

deltaP = P - nowP;
deltaQ = Q - nowQ;
% 将平衡节点和pv节点相应的量设为0 
deltaP(balance) = 0;
deltaQ(balance) = 0;
deltaQ(pvs) = 0;
unb = [deltaP;deltaQ];
maxun = max(abs(unb));
unb = sparse(unb);
end