%%__��ȡ��ǰ�����µ�Pi��Qi__
%%__����:����__
%%__�������:2019��12��15��__
%%����˵��
%% U_polar:������ĸ��ڵ��ѹ        Y:�ڵ㵼�ɾ���
% function [nowP, nowQ] = calculatepower(abs_y, nodes, angij, U, balance, pvs)
function [nowP, nowQ] = calculatepower(U_polar, Y)
% rundata ���в���
% balance, pvs ��������

%% ����
% % �������̫ռʱ���� ��취���ɾ���
% nowP = zeros(nodes,1);
% nowQ = zeros(nodes,1);
% for i = 1:nodes
%     for j = 1:nodes
%         nowP(i) = nowP(i)+U(i)*abs_y(i,j)*U(j)*cos(angij(i,j));
%         nowQ(i) = nowQ(i)+U(i)*abs_y(i,j)*U(j)*sin(angij(i,j));
%     end
% end


%% �Ľ��� 
    S = sparse(diag(U_polar)) * conj(Y * U_polar);          %�ü�����ϵ��ȡ�����ʣ����ǿα�127ҳ��4-35ʽ
    nowP = real(S);
    nowQ = imag(S);
end