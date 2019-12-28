function [sysdata, linedata, branchdata, transferdata, rundata, pvdata]=loadfile()
% 
% ���:
% sysdata: ϵͳ����
% linedata: ��·����
% jiedidata: �ӵز���
% transferdata: ��ѹ������
% rundata: ���в���
% pvdata: pv�ڵ����
% gendata: ���������
% ��һ������û��д���������
% �ĵ����ڣ�2019��12��2��

% ��ȡ�ļ�·�� fileΪ�ļ���,pathΪ·��
[file,path] = uigetfile('*.dat');
if file == 0
    return
end
% ʹ�þ���·����ȡ���� dlmread��ȡ�ָ���(�����ո��tab)
datas = dlmread([path,file]);
% ��ȡ�ָ���
index = zeros([1, 8]); % Ԥ�����ڴ�, 8��λ��
[row, ~] = size(datas);
j = 1;
% �ļ��������������0��
for i=1:row
    
    temp = datas(i,:);
    % �жϵ�ǰ���Ƿ�ȫ��
    if ~any(temp)
        index(j) = i;
        j = j + 1;
    end
end
% ����0�зָ�����
sysdata = datas(1:index(1)-1,:);            % ϵͳ����
linedata = datas(index(1)+1:index(2)-1,:);% ��·����
branchdata = datas(index(2)+1:index(3)-1,:);% �ӵز���
transferdata = datas(index(3)+1:index(4)-1,:);% ��ѹ������
rundata = datas(index(4)+1:index(5)-1,:);% ���в���
pvdata = datas(index(5)+1:index(6)-1,:);% pv�ڵ����

end