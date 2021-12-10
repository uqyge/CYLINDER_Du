function [Be,te,Ae] = calc_local_stiffness( INISTRAIN )

    global NODE_LOCAL0
    global MATE
    global PARA
    
    global GAUS1
    global WGT1
    
    global KELIST
    global INFLIST
    global INFLIST_PESUDO
    
    Be = zeros(3,6,PARA.NELEM);
    te = ones(PARA.NELEM,1);
    Ae = zeros(PARA.NELEM,1);
    
    %% Element stiffness
    for ielem=1:PARA.NELEM
        
        node1_coord=NODE_LOCAL0(ielem).X1;node2_coord=NODE_LOCAL0(ielem).X2;node3_coord=NODE_LOCAL0(ielem).X3;
        x1=node1_coord(1);y1=node1_coord(2);
        x2=node2_coord(1);y2=node2_coord(2);
        x3=node3_coord(1);y3=node3_coord(2);
        
        A2=det([1 x1 y1;
                1 x2 y2;
                1 x3 y3]); % 2A
        
        Dm = MATE.Dm;
        h = PARA.h;
        te(ielem,1) = h;
        G = MATE.EX/(1+MATE.PRXY);
        
        %% CST
         Km = zeros(9,9);
         fi = zeros(9,1);
         fi_pesudo = zeros(9,1);
         for j=1:size(GAUS1,2)
             DL = [1  0  -1;
                   0  1  -1];
             Jacobi = DL*[x1 y1;
                          x2 y2;
                          x3 y3];
             DG = Jacobi\DL;
             detJ = det(Jacobi);
             
             B = [DG(1,1)       0  0  DG(1,2)      0   0  DG(1,3)      0   0;
                  0       DG(2,1)  0  0       DG(2,2)  0  0       DG(2,3)  0;
                  DG(2,1) DG(1,1)  0  DG(2,2) DG(1,2)  0  DG(2,3) DG(1,3)  0];
              
             Km = Km + B'*Dm*B*detJ*WGT1(j)*h;
             fi = fi + B'*(Dm*INISTRAIN)*detJ*WGT1(j)*h;
             fi_pesudo = fi_pesudo + B'*(Dm*[1E-12;1E-12;0])*detJ*WGT1(j)*h;
         end
         
%          Q = [-(x2-x3)/(A2*2) -(y2-y3)/(A2*2) 1/3 -(x3-x1)/(A2*2) -(y3-y1)/(A2*2) 1/3 -(x1-x2)/(A2*2) -(y1-y2)/(A2*2) 1/3];
%          delta2 = 0.01;
%          V = A2/2*h;
%          Sr = (delta2*V*G)*(Q'*Q);
%             
%          Km = Km + Sr;
%          Km(3,3)=Km(3,3)+G*V*1E-8;
%          Km(6,6)=Km(6,6)+G*V*1E-8;
%          Km(9,9)=Km(9,9)+G*V*1E-8;
         
        %% Assemble element stiffness
         Kel = zeros(18,18);
         fel = zeros(18,1);
         fel_pesudo = zeros(18,1);
         index_b = [3 4 5 9 10 11 15 16 17];
         index_m = [1 2 6 7 8 12 13 14 18];
         
         Kel(index_b,index_b)=0;
         Kel(index_m,index_m)=Km;
         
         fel(index_b)=0;
         fel(index_m)=fi;
         
         fel_pesudo(index_b)=0;
         fel_pesudo(index_m)=fi_pesudo;
         
         KELIST(ielem).M = Kel;
         INFLIST(ielem).M = fel;
         INFLIST_PESUDO(ielem).M = fel_pesudo;
         
         B36 = [B(:,1:2),B(:,4:5),B(:,7:8)];
         
         Be(:,:,ielem) = B36;
         
         Ae(ielem,1) = 0.5*detJ;
         
    end
end

