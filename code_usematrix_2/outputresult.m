%%__�������__
%%__����:����__
%%__�������:2019��12��15��__
%% ����˵��:
%%nowP:��ǰ�й��ڵ㹦��        nowQ:��ǰ�ڵ��޹�����
%%Y:�ڵ㵼�ɾ���     alphaU:���       U_polar:�������ѹ
%%linedata:��·����
function outputresult(nowP, nowQ, U, alphaU, U_polar, Y, linedata)
    disp('���ڵ㹦��Ϊ:');
    nowS = nowP + 1i*nowQ;
    disp(nowS);
    disp('���ڵ��ѹ�����Ϊ:');
    disp([U rad2deg(alphaU)]);
    
    plot(U)
    hold on;
    plot(alphaU);
    legend('��ֵ','����');
    title('��ѹ��ֵ�����ȱ仯����ͼ');
    xlabel('�ڵ�');
    ylabel('��ֵ/����');
    % �ɿα�133ҳ 4-51ab
    
    % ���¼���Yij�Լ�Yi��Ӱ��
    S = sparse(zeros(length(Y)));
    [rows, ~] = size(linedata);
    for i = 1:rows
        linei = linedata(i, 2);
        linej = linedata(i, 3);
        S(linei, linej) = U_polar(linei)*conj((U_polar(linei) - U_polar(linej))*(-Y(linei,linej))) + U_polar(linei)*conj(U_polar(linei)*1i*linedata(i,6));
        S(linej, linei) = U_polar(linej)*conj((U_polar(linej) - U_polar(linei))*(-Y(linej,linei))) + U_polar(linej)*conj(U_polar(linej)*1i*linedata(i,6));
    end
    disp('��·����');
    disp(S);
%     hold off;
%     bar(nowP);
end