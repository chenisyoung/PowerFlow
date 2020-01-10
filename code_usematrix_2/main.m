%%__牛顿拉夫逊算法的潮流计算__
%%__作者:陈友__
%%__完成日期:2020年1月2日__
%% 变量说明
%%n:节点数        linedata:线路参数       branchdata:支路参数
%%transferdata:变压器参数
clc;clear;
%% 读取参数
[sysdata, linedata, branchdata, transferdata, rundata, pvdata]=loadfile();
    % 取出一些常用的参数
    e = sysdata(2, 1); % 迭代精度
    MaxIteration = sysdata(1, 4); % 最大迭代次数
    nodes = sysdata(1,1); % 节点数
    Sb = sysdata(1,3); %功率标幺值 基准
    balance = sysdata(3,2); % 平衡节点编号
    pvs = pvdata(:,1); % pv节点编号
tic % 开始计时
%% 形成雅可比矩阵并给相关量设初值
    Y=Y_matrix(nodes, linedata, branchdata, transferdata); % 使用复数直接计算
    [U_polar, U, alphaU] = init_variables(nodes, pvdata);  % 初始化值
%% 开始迭代求解
    loops = 0; % 记录循环次数
    while 1
        [nowP, nowQ] = calculatepower(U_polar, Y); % 先计算 当前各节点功率 即 Pi Qi
        % 求不平衡量 和 当前精度(最大不平衡量)
        [unb, maxun] = calcutedelta(nodes, rundata, nowP, nowQ, Sb, balance, pvs);
        % loss(loops+1) = maxun; % 绘制最大不平衡量下降曲线
        if maxun < e  % 如果精度足够
            disp('精度达到要求');
            disp('精度达到要求时间');
            toc % 记录计算时间
            disp(['循环次数' ': ' num2str(loops)]);
            outputresult(nowP, nowQ, U, alphaU, U_polar, Y, linedata);
            alphaU = rad2deg(alphaU); % 将弧度转化为角度
            break;
        end
        if loops > MaxIteration  % 若超过最大迭代此时则输出无法收敛
            disp('此次计算无法收敛无法收敛')
            disp(['循环次数: ' num2str(loops)]);
            break;
        end
        Jac = Jacobin(nodes, U_polar, Y, balance, pvs);  % 计算雅可比矩阵
        delta = Jac \ unb;      % 即 inv(Jac) * unb 
        alphaU = alphaU - delta(1:nodes);  % 相角的修正
        U = U - delta(nodes+1:end) .* abs(U);  % 幅值的修正
        U_polar = sparse(U.*exp(1i.*alphaU));  % 重新化为复数并稀疏
        loops = loops + 1;
    end
    disp('输出和绘图完成后时间');
    toc % 记录完成输出以及绘图之后的时间