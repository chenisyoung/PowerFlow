%%__以平启动给参数赋初始值__
%%__作者:陈友__
%%__完成日期:2019年12月15日__
%%变量说明
%% n:节点数      pvdata:pv节点数据
function [U_polar, U, alphaU] = init_variables(n, pvdata)
    U = ones(n,1);       % 初幅值1
    alphaU = zeros(n,1);  % 初相位0
    U(pvdata(:,1))=pvdata(:,2); % pv节点电压代入入 
    U_polar = U.*exp(1i.*alphaU); 
end