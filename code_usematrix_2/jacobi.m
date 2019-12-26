% �����ſɱȾ���
function jacobmatrix = jacobi(abs_y, angij, U, nodes, balance, pvs)
% abs_y ���ɷ�ֵ
% angij delta_i - delta_j - alpha ����P��QʱҲ�õ�, Ϊ�˼��ټ�����ֱ�Ӵ������
% U �����ѹ, ά��(nodes, 1)
% nodes �ڵ����
% balance ƽ��ڵ����
% pvs pv�ڵ����

%YACOBI �ü������������ſ˱Ⱦ���
U1 = diag(U);
% ʹ�þ������ Ч�ʱ�ѭ���߷ǳ���
H = U1 * (abs_y .* sin(angij) * U1)-U1 * diag(abs_y .* sin(angij) * U);                   
N = diag(abs_y .* cos(angij) * U) + U1 * (abs_y .* cos(angij)); %ʡ������ �˴�������Ui                     
J = U1 * diag(abs_y .* cos(angij) * U)-U1 * (abs_y .* cos(angij)) * U1;                    
L = diag(abs_y .* sin(angij) * U) + U1 * (abs_y .* sin(angij));                          
H(balance,:) = 0;  % ƽ��ڵ�����
H(:,balance) = 0;
H = sparse(H);
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
jacobmatrix = -[H, N;J, L]; % ����ϰ�߽�P����Q���Ϸ�, ����ͬ���ı����򼴿�

end