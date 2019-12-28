function Y=Y_matrix(sysdata, linedata, branchdata, transferdata)
   
 %% �γ���·����
    % linedata[���, �׽ڵ�, β���, ֧·����, ֧·�翹, �Եص���/2]
     Y = zeros(sysdata(1, 1));  % �������ռ���ڴ��һ��
    for t=1 : size(linedata)
        temp = linedata(t,:);
        i = temp(2);
        j = temp(3);
        com = complex(temp(4), temp(5));
        G = -1 / com; % ����
        Y(i, i) = Y(i, i) + 1i*temp(6) - G; % �Ե���
        Y(j, j) = Y(j, j) + 1i*temp(6) - G;
        Y(i, j) = G;
        Y(j, i) = G;
        if i == j
            disp('sa')
        end
    end
Y = sparse(Y);
 %% �γɽӵ�֧·����
 for i = 1:size(branchdata)
     branch_temp = branchdata(i,:);
     iline = branch_temp(1);
     Gi = 1i*branch_temp(2);
     Y(iline, iline) = Y(iline, iline) + Gi;
 end
%  Y_d = sparse(branchdata(:,1),branchdata(:,1),branchdata(:,2));
%  Y = Y + Y_d;
 %% �γɱ�ѹ������
for i = 1:size(transferdata)
    t_temp = transferdata(i,:);
    k = t_temp(6);
    Yt = 1 / complex(t_temp(4), t_temp(5));
    % i �ڵ�
    Y(t_temp(2), t_temp(2)) = Y(t_temp(2), t_temp(2)) + Yt / k^2; %#ok<*SPRIX>
    % j �ڵ�
    Y(t_temp(3), t_temp(3)) = Y(t_temp(3), t_temp(3)) + Yt;
    Y(t_temp(2), t_temp(3)) = Y(t_temp(2), t_temp(3)) - Yt / k;              % Y(t_temp(3), t_temp(3))
    Y(t_temp(3), t_temp(2)) = Y(t_temp(3), t_temp(2)) - Yt / k;
end
Y = sparse(Y);

end