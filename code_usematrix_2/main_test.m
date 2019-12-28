clc;clear;
[sysdata, linedata, branchdata, transferdata, rundata, pvdata]=loadfile();
    % 一些常用的参数
    e = sysdata(2, 1);
    maxdiedai = sysdata(1, 4);
    nodes = sysdata(1,1); % 节点数
    Sb = sysdata(1,3); %功率标幺值 基准
    balance = sysdata(3,2); % 平衡节点编号
    pvs = pvdata(:,1); % pv节点编号
    tic
    Y=Y_matrix(sysdata, linedata, branchdata, transferdata); % 使用复数直接计算 速度差不多 更容易看懂
    toc
    % Y = Form_Y_matrix(nodes, linedata(:,2:end), transferdata(:,2:end), branchdata);%这个函数使用指导书上的算法

    %% 设初值 之后弄成function
    U = ones(sysdata(1,1),1);       % 初幅值1
    alphaU = zeros(sysdata(1,1),1);  % 初相位0
    U(pvdata(:,1))=pvdata(:,2); % pv节点电压代入入 
    U_polar = U.*exp(1i.*alphaU); 
    % abs_y = abs(Y);
    %% 开始迭代求解
    toc
    loops = 0;
    while 1
        % 先计算 当前功率
        [nowP, nowQ] = calculatepower(U_polar, Y); 
        % 求不平衡量 和 当前精度
        [unb, maxun] = calcutedelta(nodes, rundata, nowP, nowQ, Sb, balance, pvs);
        if maxun < e  % 如果精度足够
            toc
            disp('精度达到要求');
            disp(['循环次数' ': ' num2str(loops)]);
            alphaU = rad2deg(alphaU);
            break;
        end
        if loops > maxdiedai
            disp('无法收敛')
            disp(['循环次数: ' num2str(loops)]);
            break;
        end
        Jac = Jacobin(nodes, U_polar, Y, balance, pvs);  % 参考他人写法
        % Jac = jacobi(abs_y, angij, U, nodes, balance, pvs); % 按公式所写 效率有点差
        % 
        delta = Jac \ unb;
        alphaU = alphaU - delta(1:nodes);  % 先是相角
        U = U - delta(nodes+1:end) .* abs(U);  % 幅值
        U_polar = U.*exp(1i.*alphaU);  % 重新化为复数
        loops = loops + 1;
    end