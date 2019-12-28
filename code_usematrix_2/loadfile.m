function [sysdata, linedata, branchdata, transferdata, rundata, pvdata]=loadfile()
% 
% 输出:
% sysdata: 系统参数
% linedata: 线路参数
% jiedidata: 接地参数
% transferdata: 变压器参数
% rundata: 运行参数
% pvdata: pv节点参数
% gendata: 发电机参数
% 进一步处理没有写在这个函数
% 文档日期：2019年12月2日

% 获取文件路径 file为文件名,path为路径
[file,path] = uigetfile('*.dat');
if file == 0
    return
end
% 使用绝对路径获取数据 dlmread读取分隔符(包括空格和tab)
datas = dlmread([path,file]);
% 获取分隔行
index = zeros([1, 8]); % 预分配内存, 8个位置
[row, ~] = size(datas);
j = 1;
% 文件最后总是有两个0行
for i=1:row
    
    temp = datas(i,:);
    % 判断当前行是否全零
    if ~any(temp)
        index(j) = i;
        j = j + 1;
    end
end
% 按照0行分隔数据
sysdata = datas(1:index(1)-1,:);            % 系统参数
linedata = datas(index(1)+1:index(2)-1,:);% 线路参数
branchdata = datas(index(2)+1:index(3)-1,:);% 接地参数
transferdata = datas(index(3)+1:index(4)-1,:);% 变压器参数
rundata = datas(index(4)+1:index(5)-1,:);% 运行参数
pvdata = datas(index(5)+1:index(6)-1,:);% pv节点参数

end