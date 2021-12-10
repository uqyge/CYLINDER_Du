function [ ] = calcPressureStiff( lambda )
    global GAUS1
    global WGT1
    global PRES
    
    global NODE2
    global ELEM
    global K
    
    eps=zeros(3,3,3);
    eps(1,2,3)=1;eps(2,3,1)=1;eps(3,1,2)=1;
    eps(1,3,2)=-1;eps(3,2,1)=-1;eps(2,1,3)=-1;
    
    Nb = [1/3 1/3 1/3];
    DLb = [1  0  -1;
           0  1  -1];
    
    for ipres=1:size(PRES,2)
        ielem=PRES(1,ipres);
        isurf=PRES(2,ipres);
        mag=-lambda*PRES(3,ipres);
        
        if(isurf==1)
            sn=[ELEM(3,ielem);ELEM(2,ielem);ELEM(1,ielem)];
        elseif(isurf==2)
            sn=[ELEM(1,ielem);ELEM(2,ielem);ELEM(3,ielem)];
        end
        
        xsurf=[NODE2(1:3,sn(1)),NODE2(1:3,sn(2)),NODE2(1:3,sn(3))];
        i0j0_surf=[sn(1)*3-2;sn(1)*3-1;sn(1)*3;sn(2)*3-2;sn(2)*3-1;sn(2)*3;sn(3)*3-2;sn(3)*3-1;sn(3)*3];
        
        elempresstiff=zeros(3*3,3*3);
        for igaus=1:size(GAUS1,2)
            dNdxi=DLb(1,:)';dNdeta=DLb(2,:)';
            dxdxi=xsurf*dNdxi;
            dxdeta=xsurf*dNdeta;
            wgt=WGT1(igaus);
            
            for inode=1:size(Nb,2)
                for jnode=1:size(Nb,2)
                    posi=3*(inode-1);posj=3*(jnode-1);
                    KK=zeros(3,3);
                    for i=1:3
                        for j=1:3
                            for k=1:3
                                item1=mag*eps(i,j,k)*Nb(inode);
                                item2=dNdxi(jnode)*dxdeta(k);
                                item3=dNdeta(jnode)*dxdxi(k);
                                KK(i,j)=KK(i,j)+item1*(item2-item3)*wgt;
                            end
                        end
                    end
                    elempresstiff([posi+1,posi+2,posi+3],[posj+1,posj+2,posj+3])=elempresstiff([posi+1,posi+2,posi+3],[posj+1,posj+2,posj+3])+KK;
                end
            end
            
        end         
        K(i0j0_surf,i0j0_surf)=K(i0j0_surf,i0j0_surf)-elempresstiff;
    end
end

                    
                    
                    
                    
    
    
    
    
    

