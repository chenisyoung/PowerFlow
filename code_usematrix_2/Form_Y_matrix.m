function Y=Form_Y_matrix(n,data_of_line,data_of_transformers,data_of_ground)

%% % % % % % % %�γɽڵ㵼�ɾ���% % % % % % % %
Gij=data_of_line(:,3)./(data_of_line(:,3).*data_of_line(:,3)+data_of_line(:,4).*data_of_line(:,4));                %��·�絼���㣺Gij=r/(r*r+x*x)                     
Bij=-data_of_line(:,4)./(data_of_line(:,3).*data_of_line(:,3)+data_of_line(:,4).*data_of_line(:,4));               %��·���ɼ��㣺Gij=r/(r*r+x*x)
% ����絼����
G=sparse(data_of_line(:,1),data_of_line(:,2),-Gij,n,n);                                    %�γɽڵ�絼����������
G=G+sparse(data_of_line(:,2),data_of_line(:,1),-Gij,n,n);                                  %�γɽڵ�絼����������
G=G+sparse(data_of_line(:,1),data_of_line(:,1),Gij,n,n);                                   %����G(i,i)
G=G+sparse(data_of_line(:,2),data_of_line(:,2),Gij,n,n);                                   %����G(j,j)��i��j�����ظ��ĵط���ͨ�����������������
% ������ɾ���
B=sparse(data_of_line(:,1),data_of_line(:,2),-Bij,n,n);                                    %�γɽڵ���ɾ���������
B=B+sparse(data_of_line(:,2),data_of_line(:,1),-Bij,n,n);                                  %�γɽڵ���ɾ���������
B=B+sparse(data_of_line(:,1),data_of_line(:,1),Bij+data_of_line(:,5),n,n);                         %����B(i,i)
B=B+sparse(data_of_line(:,2),data_of_line(:,2),Bij+data_of_line(:,5),n,n);                         %����B(j,j)��i��j�����ظ��ĵط���ͨ�����������������
% ׷�ӱ�ѹ���������ڵ����Ľ������䣬ֻ�Ǳ�ѹ���Ľ����൱����ԭ����ڵ�����һ�ӵ�֧·���ڽڵ�i��j������һ��֧·������ķ�ʽһ��   
Gij=data_of_transformers(:,3)./(data_of_transformers(:,3).*data_of_transformers(:,3)+data_of_transformers(:,4).*data_of_transformers(:,4));                  
Bij=-data_of_transformers(:,4)./(data_of_transformers(:,3).*data_of_transformers(:,3)+data_of_transformers(:,4).*data_of_transformers(:,4));
G=G+sparse(data_of_transformers(:,1),data_of_transformers(:,2),-Gij./data_of_transformers(:,5),n,n);
G=G+sparse(data_of_transformers(:,2),data_of_transformers(:,1),-Gij./data_of_transformers(:,5),n,n);
G=G+sparse(data_of_transformers(:,1),data_of_transformers(:,1),Gij./data_of_transformers(:,5)./data_of_transformers(:,5),n,n);
G=G+sparse(data_of_transformers(:,2),data_of_transformers(:,2),Gij,n,n);
B=B+sparse(data_of_transformers(:,1),data_of_transformers(:,2),-Bij./data_of_transformers(:,5),n,n);
B=B+sparse(data_of_transformers(:,2),data_of_transformers(:,1),-Bij./data_of_transformers(:,5),n,n);
B=B+sparse(data_of_transformers(:,1),data_of_transformers(:,1),Bij./data_of_transformers(:,5)./data_of_transformers(:,5),n,n);
B=B+sparse(data_of_transformers(:,2),data_of_transformers(:,2),Bij,n,n);
% ׷�ӶԵ�֧·����
B=B+sparse(data_of_ground(:,1),data_of_ground(:,1),data_of_ground(:,2),n,n);                             %���ϲ��������B��������ֱ�����
% д�ɸ�����ʽ   
Y=G+1i*B;                                                                   %�γ��Ը�����ʾ�Ľڵ㵼�ɾ���