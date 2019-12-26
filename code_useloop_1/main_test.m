clc
clear
[sysdata, linedata, branchdata, transferdata, rundata, pvdata, gendata]=loadfile();
tic
    Y=Y_matrix(sysdata, linedata, branchdata, transferdata);
    % 功率标幺值 基准
    Sb = sysdata(1,3);
    fuzhi = abs(Y);
    jiao = angle(Y);
    %% 设初值
    U = ones(sysdata(1,1),1);       % 初幅值1
    alphaU = zeros(sysdata(1,1),1);  % 初相位0
    U(pvdata(:,1))=pvdata(:,2); % pv节点电压代入入 
    %% 开始迭代求解
    %for i = 1:3 % 迭代次数
        [ddelta, dU] = Jacobian(fuzhi, jiao, sysdata, rundata, U, alphaU, pvdata, Sb);
        U = U - dU.*U;
        alphaU = alphaU - ddelta;
        % disp(['第几次迭代' ' ' num2str(i)]);
    %end
toc