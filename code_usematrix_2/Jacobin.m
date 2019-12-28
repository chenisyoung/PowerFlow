function Jb=Jacobin(n,U_polar,Y,slack_bus,PVnode)
 %% % % % % % %%形成雅克比矩阵，采用了完全极坐标条件下去计算雅克比矩阵%% % % % % % %%
        % U_polar = sparse(U_polar);写在别处 涉及矩阵运算的需都转化为稀疏矩阵
        H=sparse(imag(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))-imag(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Hij—对△Pi求△δj的偏导，Hii—对△Pi求△δi的偏导
        N=sparse(-real(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))-real(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Nij—对△Pi求△Vj的偏导，Nii—对△Pi求△Vi的偏导
        J=sparse(-real(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))+real(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Jij—对△Qi求△δj的偏导，Hii—对△Qi求△δi的偏导
        L=sparse(-imag(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))-imag(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Lij—对△Qi求△Vi的偏导，Lii—对△Qi求△Vj的偏导
        %计算H  
        H(slack_bus,:)=0;                                                          %对H中的平衡节点进行修正，将H中平衡节点对应的行置零 
        H=sparse(H);
        H(:,slack_bus)=0;                                                          %对H中的平衡节点进行修正，将H中平衡节点对应的行置零 
        H=sparse(H);
        H=H+sparse(slack_bus,slack_bus,1,n,n);                                     %H主对角元素不能为0，需要对H的对角线元素进行修改

        %计算N  
        N(:,PVnode)=0;                                                             %对N中的pv节点进行修正置零             
        N=sparse(N);
        N(slack_bus,:)=0;                                                          %对N中的平衡节点进行修正，将N中平衡节点对应的行置零                                                                 
        N=sparse(N);
        N(:,slack_bus)=0;                                                          %对N中的平衡节点进行修正，将N中平衡节点对应的列置零  
        N=sparse(N);
               
        %计算J       
        J(PVnode,:)=0;                                                            %对J中的pv节点进行修正，因为PV节点只有H和N,将L、J置零，但是J的主对角元素不能为0
        J=sparse(J);
        J(slack_bus,:)=0;                                                         %对J中的平衡节点进行修正，将J中平衡节点对应的行置零
        J=sparse(J);
        J(:,slack_bus)=0;                                                         %对J中的平衡节点进行修正，将J中平衡节点对应的列置零
        J=sparse(J);
        
        %计算L  
        L(PVnode,:)=0;                                                            %对L中的pv节点进行修正，因为PV节点只有H和N,将L、J置零，但是L的主对角元素不能为0
        L=sparse(L);
        L(:,PVnode)=0;  
        L=sparse(L);
        L=L+sparse(PVnode,PVnode,1,n,n);
        L(slack_bus,:)=0;                                                          %对L中的平衡节点进行修正，将L中平衡节点对应的行置零
        L=sparse(L);
        L(:,slack_bus)=0;                                                          %对L中的平衡节点进行修正，将L中平衡节点对应的列置零
        L=sparse(L);
        L=L+sparse(slack_bus,slack_bus,1,n,n);                                     %L主对角元素不能为0，需要对L的对角线元素进行修改
        
        %形成雅克比矩阵
        Jb=[H,N;J,L]; 