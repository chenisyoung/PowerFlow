clc
clear
[sysdata, linedata, branchdata, transferdata, rundata, pvdata, gendata]=loadfile();
    Y=Y_matrix(sysdata, linedata, branchdata, transferdata);
    % 一些常用的参数
    e = sysdata(2, 1);
    maxdiedai = sysdata(1, 4);
    nodes = sysdata(1,1); % 节点数
    Sb = sysdata(1,3); %功率标幺值 基准
    balance = sysdata(3,2); % 平衡节点编号
    pvs = pvdata(:,1); % pv节点编号
    abs_y = abs(Y); % 导纳幅值
    angle_y = angle(Y); % 导纳相角
    %% 设初值
    U = ones(sysdata(1,1),1);       % 初幅值1
    alphaU = zeros(sysdata(1,1),1);  % 初相位0
    U(pvdata(:,1))=pvdata(:,2); % pv节点电压代入入 
    %% 开始迭代求解
    loops = 1;
    while 1
        angij = getangij(angle_y, alphaU);  % 求取 i - j - y
        % 理应先计算 当前功率
        % [nowP, nowQ] = calculatepower(abs_y, nodes, alphaU, U, balance, pvs); 
        [nowP, nowQ] = calculatepower(alphaU, U, Y, balance, pvs); 
        % 求不平衡量 和 当前精度
        [unb, maxun] = calcutedelta(rundata, nowP, nowQ, Sb, balance, pvs);
        if maxun < e  % 如果精度足够
            disp(['精度达到要求' ' ' num2str(loops)]);
            break;
        end
        if loops > maxdiedai
            disp(['无法收敛' ' ' num2str(loops)]);
            break;
        end
        Jac = jacobi(abs_y, angij, U, nodes, balance, pvs);  % 求雅可比矩阵

        % 
        delta = Jac \ unb;
        alphaU = alphaU - delta(1:nodes);  % 先是相角
        U = U - delta(nodes+1:end);  % 幅值
        loops = loops + 1;
    end