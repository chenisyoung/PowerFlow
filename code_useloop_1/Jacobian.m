function [ddelta, du] = Jacobian(fuzhi, jiao, sysdata, rundata, U, alphaU, pvdata, Sb)
% U: (1, Nodes) 
%  
balance = sysdata(3,2);
Y = fuzhi;  % 下面习惯打了Y 外面也使用Y会命名混乱 故在此进行赋值操作
nodes = sysdata(1,1);
H = zeros(nodes, nodes);
N = zeros(nodes, nodes);
J = zeros(nodes, nodes);
L = zeros(nodes, nodes);

% 取出pv节点, 不计算相关J, L
pvs = pvdata(:,1);
tic % 给 循环计时
for i =1 : nodes
    
    Ht = 0;  % 临时变量 需要求和的 每个参数都有一个(ii 或 ij)
    Nt = 0;
    Lt = 0;
    Jt = 0;
    for j = 1:nodes
%         if i == 7 && j == 4
%             disp('nothing');  % 为了打个断点
%         end
        if i ~= j
            % H
            Ht = Ht + Y(i,j)*U(j)*sin(alphaU(i) - alphaU(j) - jiao(i,j)); 
            H(i,j) = -Y(i,j)*U(i)*U(j)*sin(alphaU(i) - alphaU(j) - jiao(i,j));
            % N
            N(i,j) = -Y(i,j)*U(i)*U(j)*cos(alphaU(i) - alphaU(j) - jiao(i,j));
            % J
            J(i,j) =  Y(i,j)*U(i)*U(j)*cos(alphaU(i) - alphaU(j) - jiao(i,j)); 
            Jt = Jt + Y(i,j)*U(j)*cos(alphaU(i) - alphaU(j) - jiao(i,j));
            % L 
            L(i,j) = H(i,j);  %-Y(i,j)*U(i)*U(j)*sin(alphaU(i) - alphaU(j) - jiao(i,j));
            
        end
        Nt = Nt + Y(i,j)*cos(alphaU(i) - alphaU(j) - jiao(i,j));
        Lt = Lt + Y(i,j)*U(j)*sin(alphaU(i) - alphaU(j) - jiao(i,j));
    end

    J(i,i) = -U(i)*Jt;
    L(i,i) = -(U(i)*Lt - Y(i,i)*sin(jiao(i,i))*U(i)^2);
    H(i,i) = U(i) * Ht;
    N(i,i) = -(U(i)*Nt + Y(i,i)*cos(jiao(i,i))*U(i)^2);

end
toc
    %% 拼接矩阵
    H(balance,:) = 0;  % 平衡节点置零
    H(:,balance) = 0;
    H=H+sparse(balance,balance,1,nodes,nodes);  % 对角线0修正
    H = sparse(H);
    
    N(:,pvs)=0;   
    N(balance,:) = 0;    
    N(:,balance) = 0;
    N = sparse(N);

    J(pvs,:) = 0;        %对J中的pv节点进行修正，因为PV节点只有H和N,将L、J置零，但是J的主对角元素不能为0
    J(balance,:) = 0;  %对J中的平衡节点进行修正，将J中平衡节点对应的行置零
    J(:,balance) = 0;    %对J中的平衡节点进行修正，将J中平衡节点对应的列置零
    J = sparse(J);
    
    L(pvs,:) = 0;     %对L中的pv节点进行修正，因为PV节点只有H和N,将L、J置零，但是L的主对角元素不能为0
    L(:,pvs) = 0;  
    L = L+sparse(pvs,pvs,1,nodes,nodes);
    L(balance,:) = 0;      %对L中的平衡节点进行修正，将L中平衡节点对应的行置零
    L(:,balance) = 0;      %对L中的平衡节点进行修正，将L中平衡节点对应的列置零
    L = L+sparse(balance,balance,1,nodes,nodes);     % 对角线0修正
    L = sparse(L);
        
% for i = 0:nodes-1
%         HN(:,2*i+1) = H(:,i+1);
%         HN(:,2*i+2) = N(:,i+1);
%         
%         JL(:,2*i+1) = J(:,i+1);
%         JL(:,2*i+2) = L(:,i+1);
% end
jacob = sparse([H, N;J, L]);  % 此为雅可比矩阵
%% 求出功率 及其 变化量

sumP = zeros(nodes,1);
sumQ = zeros(nodes,1);
P = rundata(:,2) - rundata(:,4); % Pi
P = P ./ Sb;
Q = rundata(:,3) - rundata(:,5); % Qi
Q = Q ./ Sb;
for i = 1:nodes
    for j = 1:nodes
        angij = alphaU(i) - alphaU(j) - jiao(i,j);
        sumP(i) = sumP(i)+U(i)*Y(i,j)*U(j)*cos(angij);
        sumQ(i) = sumQ(i)+U(i)*Y(i,j)*U(j)*sin(angij);
    end
end
deltaP = P - sumP;
deltaQ = Q - sumQ;
% 将平衡节点和pv节点相应的量设为0
deltaP(balance) = 0;
deltaQ(balance) = 0;
deltaQ(pvs) = 0;
deltasP = [deltaP;deltaQ];
%% 计算参数变化量
deltascanshu = jacob \ deltasP;

ddelta = deltascanshu(1:nodes);
du = deltascanshu(nodes+1:end);
end
