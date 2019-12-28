function Jb=Jacobin(n,U_polar,Y,slack_bus,PVnode)
 %% % % % % % %%�γ��ſ˱Ⱦ��󣬲�������ȫ������������ȥ�����ſ˱Ⱦ���%% % % % % % %%
        % U_polar = sparse(U_polar);д�ڱ� �漰����������趼ת��Ϊϡ�����
        H=sparse(imag(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))-imag(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Hij���ԡ�Pi�����j��ƫ����Hii���ԡ�Pi�����i��ƫ��
        N=sparse(-real(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))-real(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Nij���ԡ�Pi���Vj��ƫ����Nii���ԡ�Pi���Vi��ƫ��
        J=sparse(-real(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))+real(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Jij���ԡ�Qi�����j��ƫ����Hii���ԡ�Qi�����i��ƫ��
        L=sparse(-imag(sparse(diag(U_polar))*sparse(diag(conj(Y*U_polar))))-imag(sparse(diag(U_polar))*conj(Y*sparse(diag(U_polar)))));%Lij���ԡ�Qi���Vi��ƫ����Lii���ԡ�Qi���Vj��ƫ��
        %����H  
        H(slack_bus,:)=0;                                                          %��H�е�ƽ��ڵ������������H��ƽ��ڵ��Ӧ�������� 
        H=sparse(H);
        H(:,slack_bus)=0;                                                          %��H�е�ƽ��ڵ������������H��ƽ��ڵ��Ӧ�������� 
        H=sparse(H);
        H=H+sparse(slack_bus,slack_bus,1,n,n);                                     %H���Խ�Ԫ�ز���Ϊ0����Ҫ��H�ĶԽ���Ԫ�ؽ����޸�

        %����N  
        N(:,PVnode)=0;                                                             %��N�е�pv�ڵ������������             
        N=sparse(N);
        N(slack_bus,:)=0;                                                          %��N�е�ƽ��ڵ������������N��ƽ��ڵ��Ӧ��������                                                                 
        N=sparse(N);
        N(:,slack_bus)=0;                                                          %��N�е�ƽ��ڵ������������N��ƽ��ڵ��Ӧ��������  
        N=sparse(N);
               
        %����J       
        J(PVnode,:)=0;                                                            %��J�е�pv�ڵ������������ΪPV�ڵ�ֻ��H��N,��L��J���㣬����J�����Խ�Ԫ�ز���Ϊ0
        J=sparse(J);
        J(slack_bus,:)=0;                                                         %��J�е�ƽ��ڵ������������J��ƽ��ڵ��Ӧ��������
        J=sparse(J);
        J(:,slack_bus)=0;                                                         %��J�е�ƽ��ڵ������������J��ƽ��ڵ��Ӧ��������
        J=sparse(J);
        
        %����L  
        L(PVnode,:)=0;                                                            %��L�е�pv�ڵ������������ΪPV�ڵ�ֻ��H��N,��L��J���㣬����L�����Խ�Ԫ�ز���Ϊ0
        L=sparse(L);
        L(:,PVnode)=0;  
        L=sparse(L);
        L=L+sparse(PVnode,PVnode,1,n,n);
        L(slack_bus,:)=0;                                                          %��L�е�ƽ��ڵ������������L��ƽ��ڵ��Ӧ��������
        L=sparse(L);
        L(:,slack_bus)=0;                                                          %��L�е�ƽ��ڵ������������L��ƽ��ڵ��Ӧ��������
        L=sparse(L);
        L=L+sparse(slack_bus,slack_bus,1,n,n);                                     %L���Խ�Ԫ�ز���Ϊ0����Ҫ��L�ĶԽ���Ԫ�ؽ����޸�
        
        %�γ��ſ˱Ⱦ���
        Jb=[H,N;J,L]; 