function [ lambdas, strain_matrix ] = calc_stiff_intforce_Du(e1,e2,ss,DI,Ae,Be,te)
% Calculate stiffness and internal force using EICR formulation

    global KELIST
    global INFLIST
    global INFLIST_PESUDO
    global ROTMATRIX
    global ELEM
    global PARA
    global MATE

    global K
    global I
    global J
    global V
    global I0J0
    global FINT
    global BC
    global U
    global WGT1
    
    global NODE_LOCAL0
    global T0LIST
    global NODE_LOCAL
    global TLIST
    
    global RESULT

    I=zeros(PARA.NELEM*9*9,1);
    J=zeros(PARA.NELEM*9*9,1);
    V=zeros(PARA.NELEM*9*9,1);
    
    lambdas = zeros(2*PARA.NELEM,1);
    strain_matrix = zeros(2,2,PARA.NELEM);
    
    
    FINT(1:end)=0;
    
    for ielem=1:PARA.NELEM
        
        z0 = zeros(2,1);% For Lemke
        
       %% Calculate initial CR frame
        E0=T0LIST(ielem).M;
        
       %% Calculate current CR frame
        E=TLIST(ielem).M;
        
       %% Calculate local displacement and rotation
        R1g=ROTMATRIX(ELEM(1,ielem)).M;R2g=ROTMATRIX(ELEM(2,ielem)).M;R3g=ROTMATRIX(ELEM(3,ielem)).M;
        u1e=NODE_LOCAL(ielem).x1-NODE_LOCAL0(ielem).X1;u2e=NODE_LOCAL(ielem).x2-NODE_LOCAL0(ielem).X2;u3e=NODE_LOCAL(ielem).x3-NODE_LOCAL0(ielem).X3;
        R1e=E'*R1g*E0;R2e=E'*R2g*E0;R3e=E'*R3g*E0;
        Omega1e=calc_Omega(R1e);Omega2e=calc_Omega(R2e);Omega3e=calc_Omega(R3e);
        theta1e=[Omega1e(3,2);Omega1e(1,3);Omega1e(2,1)];theta2e=[Omega2e(3,2);Omega2e(1,3);Omega2e(2,1)];theta3e=[Omega3e(3,2);Omega3e(1,3);Omega3e(2,1)];
        pe=[u1e;theta1e;u2e;theta2e;u3e;theta3e];
        
%        %% Complementary problem
%         pe_mem = [u1e(1:2);u2e(1:2);u3e(1:2)];
%         strain = Be(:,:,ielem)*pe_mem;
%         strain_matrix(:,:,ielem) = [strain(1,1),strain(3,1)/2;...
%                          strain(3,1)/2,strain(2,1)];
%         [ VV, DD ] = eig(strain_matrix(:,:,ielem));
%         l1=VV(1,1);m1=VV(2,1);
%         l2=VV(1,2);m2=VV(2,2);
%         strain_p = [DD(1,1);DD(2,2)];
%         Q = [l1^2,m1^2,l1*m1;l2^2,m2^2,l2*m2];
%         He = ss*DI*Q*Be(:,:,ielem);
%         We = ss*te(ielem)*Ae(ielem)*Be(:,:,ielem)'*Q'*DI;
%         qqe = -He*pe_mem;
%         AAe = ss*(e1*e2/(e2-e1))*[1,0;0,1] + DI;
%         [ lambdae ] =lemke(AAe,qqe,z0);
%         lambdas(2*ielem-1:2*ielem,1) = lambdae;
%         Tmp1 = We*lambdae;
%         Tmp = [Tmp1(1:2);0;0;0;0;Tmp1(3:4);0;0;0;0;Tmp1(5:6);0;0;0;0];

        %% Zonglinag-Du's Method
        % Calculate Principle Strain
        pe_mem = [u1e(1:2);u2e(1:2);u3e(1:2)];
        strain = Be(:,:,ielem)*pe_mem;
        strain_matrix(:,:,ielem) =  [strain(1,1),strain(3,1)/2;
                                     strain(3,1)/2,strain(2,1)];
        [VV ,DD] = eig(strain_matrix(:,:,ielem));
        l1=VV(1,1);m1=VV(2,1);
        l2=VV(1,2);m2=VV(2,2);
        strain_p = [DD(1,1);DD(2,2)];
        % Calculate Subspace D
        E_m = e1;
        E_p = e2;
        mu_m = MATE.PRXY;
        mu_p = mu_m*E_p/E_m;
        
%         D1p = [-E_p/(C0^2*E_p^2-1) -C0*E_p^2/(C0^2*E_p^2-1);
%                -C0*E_p^2/(C0^2*E_p^2-1) -E_p/(C0^2*E_p^2-1)];
%         D2p = [-E_m/(C0^2*E_m^2-1) -C0*E_m^2/(C0^2*E_m^2-1);
%                -C0*E_m^2/(C0^2*E_m^2-1) -E_m/(C0^2*E_m^2-1)];
%         D3p = [-E_p/(C0^2*E_m*E_p-1) -C0*E_p*E_m/(C0^2*E_m*E_p-1);
%                -C0*E_p*E_m/(C0^2*E_m*E_p-1) -E_m/(C0^2*E_m*E_p-1)];
%         D4p = [-E_m/(C0^2*E_m*E_p-1) -C0*E_p*E_m/(C0^2*E_m*E_p-1);
%                -C0*E_p*E_m/(C0^2*E_m*E_p-1) -E_p/(C0^2*E_m*E_p-1)];

        D1p = 1/(1-mu_m*mu_m)*[E_m,E_m*mu_m;E_m*mu_m,E_m];
        D2p = 1/(1-mu_p*mu_p)*[E_p,E_p*mu_p;E_p*mu_p,E_p];
        D3p = 1/(1-mu_m*mu_p)*[E_m,E_m*mu_p;E_p*mu_m,E_p];
        D4p = 1/(1-mu_m*mu_p)*[E_p,E_m*mu_p;E_p*mu_m,E_m];

%         D1p = 1/((1+mu_m)^2)*[E_m,E_m*mu_m;E_m*mu_m,E_m];
%         D2p = 1/((1+mu_p)^2)*[E_p,E_p*mu_p;E_p*mu_p,E_p];
%         D3p = 1/(1-mu_m*mu_p)*[E_m,E_m*mu_p;E_m*mu_p,E_p];
%         D4p = 1/(1-mu_m*mu_p)*[E_p,E_m*mu_p;E_m*mu_p,E_m];


        stress = zeros(3,1);
        Dtan = zeros(3,3);
        % Calculate Auxiliary Quantities
        n1 = [l1;m1];
        n2 = [l2;m2];
         
        n11 = VV(1,1); n12 = VV(2,1); n21 = VV(1,2); n22 = VV(2,2);
        Q = [n11^2,n12^2,n11*n12; n21^2,n22^2,n21*n22];
        e_e1 = [1,0;0,0]; e_e2 = [0,0;0,1]; e_e3 = 0.5*[0,1;1,0];
        P = strain_p;
        
        n1_e1 = -1/(P(1)-P(2))*n2'*((n1'*e_e1*n1)*n1-e_e1*n1)*n2;
        n2_e1 = 1/(P(1)-P(2))*n1'*((n2'*e_e1*n2)*n2-e_e1*n2)*n1;
        n1_e2 = -1/(P(1)-P(2))*n2'*((n1'*e_e2*n1)*n1-e_e2*n1)*n2;
        n2_e2 = 1/(P(1)-P(2))*n1'*((n2'*e_e2*n2)*n2-e_e2*n2)*n1;
        n1_e3 = -1/(P(1)-P(2))*n2'*((n1'*e_e3*n1)*n1-e_e3*n1)*n2;
        n2_e3 = 1/(P(1)-P(2))*n1'*((n2'*e_e3*n2)*n2-e_e3*n2)*n1;
       
        Q_e1 = [2*n11*n1_e1(1),2*n12*n1_e1(2),n11*n1_e1(2)+n12*n1_e1(1);
                2*n21*n2_e1(1),2*n22*n2_e1(2),n21*n2_e1(2)+n22*n2_e1(1)];
        Q_e2 = [2*n11*n1_e2(1),2*n12*n1_e2(2),n11*n1_e2(2)+n12*n1_e2(1);
                2*n21*n2_e2(1),2*n22*n2_e2(2),n21*n2_e2(2)+n22*n2_e2(1)];
        Q_e3 = [2*n11*n1_e3(1),2*n12*n1_e3(2),n11*n1_e3(2)+n12*n1_e3(1);
                2*n21*n2_e3(1),2*n22*n2_e3(2),n21*n2_e3(2)+n22*n2_e3(1)];
            
%         % Check Subspace
%         if(strain_p(1)==strain_p(2) && (strain_p(1)+mu_p*strain_p(2))>0 && (strain_p(2)+mu_p*strain_p(1))>0)
        if(strain_p(1)==strain_p(2) && (strain_p(1)+mu_m*strain_p(2))>0 && (strain_p(2)+mu_m*strain_p(1))>0)
            Dtan = MATE.Dm;
            stress = Dtan*strain;
%         elseif(strain_p(1)==strain_p(2) && (strain_p(1)+mu_m*strain_p(2))<=0 && (strain_p(2)+mu_m*strain_p(1))<=0)
        elseif(strain_p(1)==strain_p(2) && (strain_p(1)+mu_p*strain_p(2))<=0 && (strain_p(2)+mu_p*strain_p(1))<=0)
            Dtan = E_p/(1-mu_p^2)*[1 mu_p 0;
                                   mu_p 1 0;
                                   0  0 (1-mu_p)/2];
            stress = Dtan*strain;
%         elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_p*strain_p(2))>0 && (strain_p(2)+mu_p*strain_p(1))>0) %S1 Subspace
        elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_m*strain_p(2))>0 && (strain_p(2)+mu_m*strain_p(1))>0) %S1 Subspace
            D1 = Q'*D1p*Q;
            Dtan(1,1) = D1(1,1) + Q_e1(:,1)'*D1p*Q*strain + Q(:,1)'*D1p*Q_e1*strain;
            Dtan(1,2) = D1(1,2) + Q_e2(:,1)'*D1p*Q*strain + Q(:,1)'*D1p*Q_e2*strain;
            Dtan(1,3) = D1(1,3) + Q_e3(:,1)'*D1p*Q*strain + Q(:,1)'*D1p*Q_e3*strain;
            Dtan(2,1) = D1(2,1) + Q_e1(:,2)'*D1p*Q*strain + Q(:,2)'*D1p*Q_e1*strain;
            Dtan(2,2) = D1(2,2) + Q_e2(:,2)'*D1p*Q*strain + Q(:,2)'*D1p*Q_e2*strain;
            Dtan(2,3) = D1(2,3) + Q_e3(:,2)'*D1p*Q*strain + Q(:,2)'*D1p*Q_e3*strain;
            Dtan(3,1) = D1(3,1) + Q_e1(:,3)'*D1p*Q*strain + Q(:,3)'*D1p*Q_e1*strain;
            Dtan(3,2) = D1(3,2) + Q_e2(:,3)'*D1p*Q*strain + Q(:,3)'*D1p*Q_e2*strain;
            Dtan(3,3) = D1(3,3) + Q_e3(:,3)'*D1p*Q*strain + Q(:,3)'*D1p*Q_e3*strain;           
            stress = D1*strain;
%         elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_m*strain_p(2))<=0 && (strain_p(2)+mu_m*strain_p(1))<=0) %S2 Subspace
        elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_p*strain_p(2))<=0 && (strain_p(2)+mu_p*strain_p(1))<=0) %S2 Subspace
            D2 = Q'*D2p*Q;
            Dtan(1,1) = D2(1,1) + Q_e1(:,1)'*D2p*Q*strain + Q(:,1)'*D2p*Q_e1*strain;
            Dtan(1,2) = D2(1,2) + Q_e2(:,1)'*D2p*Q*strain + Q(:,1)'*D2p*Q_e2*strain;
            Dtan(1,3) = D2(1,3) + Q_e3(:,1)'*D2p*Q*strain + Q(:,1)'*D2p*Q_e3*strain;
            Dtan(2,1) = D2(2,1) + Q_e1(:,2)'*D2p*Q*strain + Q(:,2)'*D2p*Q_e1*strain;
            Dtan(2,2) = D2(2,2) + Q_e2(:,2)'*D2p*Q*strain + Q(:,2)'*D2p*Q_e2*strain;
            Dtan(2,3) = D2(2,3) + Q_e3(:,2)'*D2p*Q*strain + Q(:,2)'*D2p*Q_e3*strain;
            Dtan(3,1) = D2(3,1) + Q_e1(:,3)'*D2p*Q*strain + Q(:,3)'*D2p*Q_e1*strain;
            Dtan(3,2) = D2(3,2) + Q_e2(:,3)'*D2p*Q*strain + Q(:,3)'*D2p*Q_e2*strain;
            Dtan(3,3) = D2(3,3) + Q_e3(:,3)'*D2p*Q*strain + Q(:,3)'*D2p*Q_e3*strain;           
            stress = D2*strain;
%         elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_m*strain_p(2))>0 && (strain_p(2)+mu_p*strain_p(1))<=0) %S3 Subspace
        elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_p*strain_p(2))>0 && (strain_p(2)+mu_m*strain_p(1))<=0) %S3 Subspace
            D3 = Q'*D3p*Q;
            Dtan(1,1) = D3(1,1) + Q_e1(:,1)'*D3p*Q*strain + Q(:,1)'*D3p*Q_e1*strain;
            Dtan(1,2) = D3(1,2) + Q_e2(:,1)'*D3p*Q*strain + Q(:,1)'*D3p*Q_e2*strain;
            Dtan(1,3) = D3(1,3) + Q_e3(:,1)'*D3p*Q*strain + Q(:,1)'*D3p*Q_e3*strain;
            Dtan(2,1) = D3(2,1) + Q_e1(:,2)'*D3p*Q*strain + Q(:,2)'*D3p*Q_e1*strain;
            Dtan(2,2) = D3(2,2) + Q_e2(:,2)'*D3p*Q*strain + Q(:,2)'*D3p*Q_e2*strain;
            Dtan(2,3) = D3(2,3) + Q_e3(:,2)'*D3p*Q*strain + Q(:,2)'*D3p*Q_e3*strain;
            Dtan(3,1) = D3(3,1) + Q_e1(:,3)'*D3p*Q*strain + Q(:,3)'*D3p*Q_e1*strain;
            Dtan(3,2) = D3(3,2) + Q_e2(:,3)'*D3p*Q*strain + Q(:,3)'*D3p*Q_e2*strain;
            Dtan(3,3) = D3(3,3) + Q_e3(:,3)'*D3p*Q*strain + Q(:,3)'*D3p*Q_e3*strain;           
            stress = D3*strain;
%         elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_p*strain_p(2))<=0 && (strain_p(2)+mu_m*strain_p(1))>0) %S4 Subspace
        elseif(strain_p(1)~=strain_p(2) && (strain_p(1)+mu_m*strain_p(2))<=0 && (strain_p(2)+mu_p*strain_p(1))>0) %S4 Subspace
            D4 = Q'*D4p*Q;
            Dtan(1,1) = D4(1,1) + Q_e1(:,1)'*D4p*Q*strain + Q(:,1)'*D4p*Q_e1*strain;
            Dtan(1,2) = D4(1,2) + Q_e2(:,1)'*D4p*Q*strain + Q(:,1)'*D4p*Q_e2*strain;
            Dtan(1,3) = D4(1,3) + Q_e3(:,1)'*D4p*Q*strain + Q(:,1)'*D4p*Q_e3*strain;
            Dtan(2,1) = D4(2,1) + Q_e1(:,2)'*D4p*Q*strain + Q(:,2)'*D4p*Q_e1*strain;
            Dtan(2,2) = D4(2,2) + Q_e2(:,2)'*D4p*Q*strain + Q(:,2)'*D4p*Q_e2*strain;
            Dtan(2,3) = D4(2,3) + Q_e3(:,2)'*D4p*Q*strain + Q(:,2)'*D4p*Q_e3*strain;
            Dtan(3,1) = D4(3,1) + Q_e1(:,3)'*D4p*Q*strain + Q(:,3)'*D4p*Q_e1*strain;
            Dtan(3,2) = D4(3,2) + Q_e2(:,3)'*D4p*Q*strain + Q(:,3)'*D4p*Q_e2*strain;
            Dtan(3,3) = D4(3,3) + Q_e3(:,3)'*D4p*Q*strain + Q(:,3)'*D4p*Q_e3*strain;           
            stress = D4*strain;
        end
        %% 计算主应力，判断是否褶皱
        stm = [stress(1) stress(3);
               stress(3) stress(2)];
        [~,ds] = eig(stm);
        if(ds(1,1)*ds(2,2)<0)
            RESULT.ISWRINKLE = 1;
        end
        % Calculate Internal Force and Tangent Stiffness
        node1_coord = NODE_LOCAL0(ielem).X1;node2_coord = NODE_LOCAL0(ielem).X2;node3_coord = NODE_LOCAL0(ielem).X3;
        x1=node1_coord(1);y1=node1_coord(2);
        x2=node2_coord(1);y2=node2_coord(2);
        x3=node3_coord(1);y3=node3_coord(2);
        h = PARA.h;
        DL = [1  0  -1;
              0  1  -1];
        Jacobi = DL*[x1 y1;
                     x2 y2;
                     x3 y3];
        detJ = det(Jacobi);
        Ke = zeros(18,18);
        fe = zeros(18,1);
        index_m = [1 2 7 8 13 14];
        Ke(index_m,index_m) = Be(:,:,ielem)'*Dtan*Be(:,:,ielem)*detJ*WGT1(1)*h;  
        fe(index_m) = Be(:,:,ielem)'*stress*detJ*WGT1(1)*h;
        
       %% Calculate H matrix
        H = zeros(18,18);
        eta1 = calc_eta(theta1e);eta2 = calc_eta(theta2e);eta3 = calc_eta(theta3e);
        H(1:3,1:3)=eye(3);H(7:9,7:9)=eye(3);H(13:15,13:15)=eye(3);
        H(4:6,4:6)=eye(3)-0.5*Omega1e+eta1*Omega1e*Omega1e;
        H(10:12,10:12)=eye(3)-0.5*Omega2e+eta2*Omega2e*Omega2e;
        H(16:18,16:18)=eye(3)-0.5*Omega3e+eta3*Omega3e*Omega3e;
        
       %% Calcualte projection matrix (S and G matrix)
        x2e=NODE_LOCAL(ielem).x2(1);x3e=NODE_LOCAL(ielem).x3(1);y3e=NODE_LOCAL(ielem).x3(2);
        S13=[0 0    0;
             0 0 -x2e;
             0 x2e 0];
        S15=[0     0   y3e;
             0     0  -x3e;
             -y3e x3e    0];
        G11=[0 0 (x3e-x2e)/(y3e*x2e);
             0 0               1/x2e;
             0 -1/x2e              0];
        G13=[0 0 -x3e/(y3e*x2e);
             0 0 -1/x2e;
             0 1/x2e 0];
        G15=[0 0 1/y3e;
             0 0 0;
             0 0 0];
        
        S=[zeros(3,3) eye(3) S13 eye(3) S15 eye(3)]';
        G=[G11 zeros(3,3) G13 zeros(3,3) G15 zeros(3,3)];
        P = eye(18,18)-S*G;
        
       %% Calculate Corotation matrix T
        T = H*P;
        Ediag = zeros(18,18);
        Ediag(1:3,1:3)=E;Ediag(4:6,4:6)=E;Ediag(7:9,7:9)=E;Ediag(10:12,10:12)=E;Ediag(13:15,13:15)=E;Ediag(16:18,16:18)=E;
        
       %% Local element calculation
        fe_pesudo = fe + INFLIST_PESUDO(ielem).M; % Pesudo Internal Force
        
       %% Global internal force
        fg = Ediag*T'*fe;
        
       %% Global material stiffness
        Km = T'*Ke*T;
        
       %% Global stress stiffness
        fnm = T'*fe_pesudo;
        Fnm = [form_S(fnm(1:3));
               form_S(fnm(4:6));
               form_S(fnm(7:9));
               form_S(fnm(10:12));
               form_S(fnm(13:15));
               form_S(fnm(16:18))];
        Kgr = -Fnm*G;
        
        Fn = -[form_S(fnm(1:3))';
               zeros(3,3);
               form_S(fnm(7:9))';
               zeros(3,3);
               form_S(fnm(13:15))';
               zeros(3,3)];
        Kgp = -G'*Fn'*P;
        
       %% Tangent stiffness in global system
        Kg = Ediag*(Km+Kgr+Kgp)*Ediag';

       %% Membrane tangent stiffness in global system
        memdof = [1;2;3;7;8;9;13;14;15];
        Kgmem = Kg(memdof,memdof);
        
       %% Membrane internal force
        fgmem = fg(memdof);

       %% Assemble
        for m=1:9
            for n=1:9
                posimn=(ielem-1)*81+(m-1)*9+n;
                I(posimn)=I0J0(m,ielem);
                J(posimn)=I0J0(n,ielem);
                V(posimn)=Kgmem(m,n);
            end
        end
        FINT(I0J0(:,ielem))=FINT(I0J0(:,ielem))+fgmem;
    end
    
    K = sparse(I,J,V,PARA.NNODE*3,PARA.NNODE*3);
    
    % Penalty for BC
    for i=1:PARA.NBC
        posi=BC(1,i);value=BC(2,i);
        FINT(posi)=FINT(posi)+1E8*K(posi,posi)*(U(posi)-value);   
    end
    
    for i=1:PARA.NBC
        posi=BC(1,i);
        K(posi,posi)=K(posi,posi)*(1+1E8);
    end
end

