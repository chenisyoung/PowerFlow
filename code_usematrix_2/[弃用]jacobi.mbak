function jacobmatrix = jacobi(abs_y, angij, U, nodes, balance, pvs)
% 返回雅可比矩阵
% abs_y 导纳幅值
% angij delta_i - delta_j - alpha 
% U 各点电压, 维度(nodes, 1)
% nodes 节点个数
% balance 平衡节点序号
% pvs pv节点序号

%用极坐标来计算雅克比矩阵
U1 = diag(U);
% 使用矩阵相乘 效率比循环高非常多
H = U1 * (abs_y .* sin(angij) * U1)-U1 * diag(abs_y .* sin(angij) * U);                   
N = diag(abs_y .* cos(angij) * U) + U1 * (abs_y .* cos(angij)); %省掉计算 此处不乘以Ui                     
J = U1 * diag(abs_y .* cos(angij) * U)-U1 * (abs_y .* cos(angij)) * U1;                    
L = diag(abs_y .* sin(angij) * U) + U1 * (abs_y .* sin(angij));   

H(balance,:) = 0;  % 平衡节点置零
H(:,balance) = 0;
H = sparse(H);

N(:,pvs)=0;   
N(balance,:) = 0;    
N(:,balance) = 0;
N = sparse(N);

J(pvs,:) = 0;        %对J中的pv节点进行修正，因为PV节点只有H和N,将L、J置零，但是J的主对角元素不能为0
J(balance,:) = 0;  %对J中的平衡节点进行修正，将J中平衡节点对应的行置零
J(:,balance) = 0;    %对J中的平衡节点进行修正，将J中平衡节点对应的列置零
J = sparse(J);

L(pvs,:) = 0;     %对L中的pv节点进行修正，因为PV节点只有H和N,将L、J置零
L(:,pvs) = 0;  
% L = L+sparse(pvs,pvs,1,nodes,nodes);
L(balance,:) = 0;      %对L中的平衡节点进行修正，将L中平衡节点对应的行置零
L(:,balance) = 0;      %对L中的平衡节点进行修正，将L中平衡节点对应的列置零
L = sparse(L);
jacobmatrix = [H, N;J, L]; % 个人习惯将P排在Q的上方, 后面同样改变排序即可
        a = zeros(nodes*2);
        jacobmatrix=sparse(jacobmatrix);
        for i = 1:nodes*2
            if jacobmatrix(i,i) == 0
                a(i,i) = 1;
            end
        end
        jacobmatrix = jacobmatrix + a;

end