
function Y=Y_matrix(sysdata, linedata, branchdata, transferdata)
   
 %% 形成线路参数
    % linedata[编号, 首节点, 尾结点, 支路电阻, 支路电抗, 对地导纳/2]
     Y = zeros(sysdata(1, 1));  % 试验后发现这个比下面的快 但是可能占用内存多一点
    for t=1 : size(linedata)
        temp = linedata(t,:);
        i = temp(2);
        j = temp(3);
        com = complex(temp(4), temp(5));
        % com = temp(4) + 1i*temp(5);
        G = -1 / com; % 导纳
        Y(i, i) = Y(i, i) + 1i*temp(6) - G; % 自导纳
        Y(j, j) = Y(j, j) + 1i*temp(6) - G;
        Y(i, j) = G;
        Y(j, i) = G;
    end
    % Y = sparse(Y);
    % 使用矩阵
%     G = linedata(:, 4);
%     B = linedata(:, 5);
%     hdaona = -1 ./ complex(G, B);  % 使用点除 得到互导纳
%     Y = sparse(linedata(:,2),linedata(:,3),hdaona);
%     for k=1 : size(linedata)
%         Y(linedata(k, 2),linedata(k, 2)) = Y(linedata(k, 2),linedata(k, 2)) + linedata(k, 6);
%         Y(linedata(k, 3),linedata(k, 3)) = Y(linedata(k, 3),linedata(k, 3)) + linedata(k, 6);
%     end
 %% 形成接地支路参数
 for i = 1:size(branchdata)
     branch_temp = branchdata(i,:);
     iline = branch_temp(1);
     Gi = 1i*branch_temp(2);
     Y(iline, iline) = Y(iline, iline) + Gi;
 end
%  Y_d = sparse(branchdata(:,1),branchdata(:,1),branchdata(:,2));
%  Y = Y + Y_d;
 %% 形成变压器参数
for i = 1:size(transferdata)
    t_temp = transferdata(i,:);
    k = t_temp(6);
    Yt = 1 / complex(t_temp(4), t_temp(5));
    % i 节点
    Y(t_temp(2), t_temp(2)) = Y(t_temp(2), t_temp(2)) + Yt / k^2; %#ok<*SPRIX>
    % j 节点
    Y(t_temp(3), t_temp(3)) = Y(t_temp(3), t_temp(3)) + Yt;
    Y(t_temp(2), t_temp(3)) = Y(t_temp(2), t_temp(3)) - Yt / k;              % Y(t_temp(3), t_temp(3))
    Y(t_temp(3), t_temp(2)) = Y(t_temp(3), t_temp(2)) - Yt / k;
end
Y = sparse(Y);

end