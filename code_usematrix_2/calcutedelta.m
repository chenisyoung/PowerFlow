%%__计算不平衡量__
%%__作者:陈友__
%%__完成日期:2019年12月15日__
%%变量说明
%% n:节点数        rundata:运行参数,即各节点功率     
%% nowP:当前参数计算的各节点有功功率
%% nowQ:当前参数计算的各节点无功功率
%% Sb:基准容量      balance:平衡节点编号      pvs:平衡节点编号
%% unb:各节点不平衡量      maxun:最大不平衡量
function [unb, maxun] = calcutedelta(n, rundata, nowP, nowQ, Sb, balance, pvs)
% 这里有个要注意的地方 14节点和30节点的PQ都是按节点序号给的 400的和1047的打乱了
% 之前没有取序号出来 直接默认按序号,故出错 找了很久bug才找到
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