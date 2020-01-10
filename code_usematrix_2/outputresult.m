%%__输出数据__
%%__作者:陈友__
%%__完成日期:2019年12月15日__
%% 变量说明:
%%nowP:当前有功节点功率        nowQ:当前节点无功功率
%%Y:节点导纳矩阵     alphaU:相角       U_polar:极坐标电压
%%linedata:线路参数
function outputresult(nowP, nowQ, U, alphaU, U_polar, Y, linedata)
    disp('各节点功率为:');
    nowS = nowP + 1i*nowQ;
    disp(nowS);
    disp('各节点电压及相角为:');
    disp([U rad2deg(alphaU)]);
    
    plot(U)
    hold on;
    plot(alphaU);
    legend('幅值','弧度');
    title('电压幅值及弧度变化曲线图');
    xlabel('节点');
    ylabel('幅值/弧度');
    % 由课本133页 4-51ab
    
    % 以下计算Yij以及Yi的影响
    S = sparse(zeros(length(Y)));
    [rows, ~] = size(linedata);
    for i = 1:rows
        linei = linedata(i, 2);
        linej = linedata(i, 3);
        S(linei, linej) = U_polar(linei)*conj((U_polar(linei) - U_polar(linej))*(-Y(linei,linej))) + U_polar(linei)*conj(U_polar(linei)*1i*linedata(i,6));
        S(linej, linei) = U_polar(linej)*conj((U_polar(linej) - U_polar(linei))*(-Y(linej,linei))) + U_polar(linej)*conj(U_polar(linej)*1i*linedata(i,6));
    end
    disp('线路功率');
    disp(S);
%     hold off;
%     bar(nowP);
end