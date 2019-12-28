% function [nowP, nowQ] = calculatepower(abs_y, nodes, angij, U, balance, pvs)
function [nowP, nowQ] = calculatepower(U_polar, Y, balance, pvs)
% rundata 运行参数
% balance, pvs 用来置零

%% 弃用
% % 这个函数太占时间了 想办法换成矩阵
% nowP = zeros(nodes,1);
% nowQ = zeros(nodes,1);
% for i = 1:nodes
%     for j = 1:nodes
%         nowP(i) = nowP(i)+U(i)*abs_y(i,j)*U(j)*cos(angij(i,j));
%         nowQ(i) = nowQ(i)+U(i)*abs_y(i,j)*U(j)*sin(angij(i,j));
%     end
% end


%% 改进后 
    S = sparse(diag(U_polar)) * conj(Y * U_polar);          %用极坐标系求取复功率，就是课本127页的4-35式
    nowP = real(S);
    nowQ = imag(S);
    
%     nowP(balance) = 0;
% 	nowQ(balance) = 0;
% 	nowQ(pvs) = 0;
end