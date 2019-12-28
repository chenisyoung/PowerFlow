function Y=Form_Y_matrix(n,data_of_line,data_of_transformers,data_of_ground)

%% % % % % % % %形成节点导纳矩阵% % % % % % % %
Gij=data_of_line(:,3)./(data_of_line(:,3).*data_of_line(:,3)+data_of_line(:,4).*data_of_line(:,4));                %线路电导计算：Gij=r/(r*r+x*x)                     
Bij=-data_of_line(:,4)./(data_of_line(:,3).*data_of_line(:,3)+data_of_line(:,4).*data_of_line(:,4));               %线路电纳计算：Gij=r/(r*r+x*x)
% 计算电导矩阵
G=sparse(data_of_line(:,1),data_of_line(:,2),-Gij,n,n);                                    %形成节点电导矩阵上三角
G=G+sparse(data_of_line(:,2),data_of_line(:,1),-Gij,n,n);                                  %形成节点电导矩阵下三角
G=G+sparse(data_of_line(:,1),data_of_line(:,1),Gij,n,n);                                   %计算G(i,i)
G=G+sparse(data_of_line(:,2),data_of_line(:,2),Gij,n,n);                                   %计算G(j,j)（i和j会有重复的地方，通过这两步可以算出）
% 计算电纳矩阵
B=sparse(data_of_line(:,1),data_of_line(:,2),-Bij,n,n);                                    %形成节点电纳矩阵上三角
B=B+sparse(data_of_line(:,2),data_of_line(:,1),-Bij,n,n);                                  %形成节点电纳矩阵下三角
B=B+sparse(data_of_line(:,1),data_of_line(:,1),Bij+data_of_line(:,5),n,n);                         %计算B(i,i)
B=B+sparse(data_of_line(:,2),data_of_line(:,2),Bij+data_of_line(:,5),n,n);                         %计算B(j,j)（i和j会有重复的地方，通过这两步可以算出）
% 追加变压器参数，节点矩阵的阶数不变，只是变压器的接入相当于在原网络节点增加一接地支路，在节点i，j增加了一条支路与上面的方式一样   
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
% 追加对地支路参数
B=B+sparse(data_of_ground(:,1),data_of_ground(:,1),data_of_ground(:,2),n,n);                             %在上步中已算出B，在这里直接相加
% 写成复数形式   
Y=G+1i*B;                                                                   %形成以复数表示的节点导纳矩阵