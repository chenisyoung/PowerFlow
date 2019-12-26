function [ddelta, du] = Jacobian(fuzhi, jiao, sysdata, rundata, U, alphaU, pvdata, Sb)
% U: (1, Nodes) 
%  
balance = sysdata(3,2);
Y = fuzhi;  % ����ϰ�ߴ���Y ����Ҳʹ��Y���������� ���ڴ˽��и�ֵ����
nodes = sysdata(1,1);
H = zeros(nodes, nodes);
N = zeros(nodes, nodes);
J = zeros(nodes, nodes);
L = zeros(nodes, nodes);

% ȡ��pv�ڵ�, ���������J, L
pvs = pvdata(:,1);
tic % �� ѭ����ʱ
for i =1 : nodes
    
    Ht = 0;  % ��ʱ���� ��Ҫ��͵� ÿ����������һ��(ii �� ij)
    Nt = 0;
    Lt = 0;
    Jt = 0;
    for j = 1:nodes
%         if i == 7 && j == 4
%             disp('nothing');  % Ϊ�˴���ϵ�
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
    %% ƴ�Ӿ���
    H(balance,:) = 0;  % ƽ��ڵ�����
    H(:,balance) = 0;
    H=H+sparse(balance,balance,1,nodes,nodes);  % �Խ���0����
    H = sparse(H);
    
    N(:,pvs)=0;   
    N(balance,:) = 0;    
    N(:,balance) = 0;
    N = sparse(N);

    J(pvs,:) = 0;        %��J�е�pv�ڵ������������ΪPV�ڵ�ֻ��H��N,��L��J���㣬����J�����Խ�Ԫ�ز���Ϊ0
    J(balance,:) = 0;  %��J�е�ƽ��ڵ������������J��ƽ��ڵ��Ӧ��������
    J(:,balance) = 0;    %��J�е�ƽ��ڵ������������J��ƽ��ڵ��Ӧ��������
    J = sparse(J);
    
    L(pvs,:) = 0;     %��L�е�pv�ڵ������������ΪPV�ڵ�ֻ��H��N,��L��J���㣬����L�����Խ�Ԫ�ز���Ϊ0
    L(:,pvs) = 0;  
    L = L+sparse(pvs,pvs,1,nodes,nodes);
    L(balance,:) = 0;      %��L�е�ƽ��ڵ������������L��ƽ��ڵ��Ӧ��������
    L(:,balance) = 0;      %��L�е�ƽ��ڵ������������L��ƽ��ڵ��Ӧ��������
    L = L+sparse(balance,balance,1,nodes,nodes);     % �Խ���0����
    L = sparse(L);
        
% for i = 0:nodes-1
%         HN(:,2*i+1) = H(:,i+1);
%         HN(:,2*i+2) = N(:,i+1);
%         
%         JL(:,2*i+1) = J(:,i+1);
%         JL(:,2*i+2) = L(:,i+1);
% end
jacob = sparse([H, N;J, L]);  % ��Ϊ�ſɱȾ���
%% ������� ���� �仯��

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
% ��ƽ��ڵ��pv�ڵ���Ӧ������Ϊ0
deltaP(balance) = 0;
deltaQ(balance) = 0;
deltaQ(pvs) = 0;
deltasP = [deltaP;deltaQ];
%% ��������仯��
deltascanshu = jacob \ deltasP;

ddelta = deltascanshu(1:nodes);
du = deltascanshu(nodes+1:end);
end
