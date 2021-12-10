function [ Fpres ] = calcPressureForce( lambda )

    global GAUS1
    global WGT1
    global PRES
    global PARA
    
    global NODE2
    global ELEM
    
    Nb = [1/3 1/3 1/3];
    DLb = [1  0  -1;
           0  1  -1];
       
    Fpres = zeros(PARA.NNODE*3,1);
    
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
        
        elempres=zeros(3*3,1);
        for igaus=1:size(GAUS1,2)
            dNdxi=DLb(1,:)';dNdeta=DLb(2,:)';
            dxdxi=xsurf*dNdxi;
            dxdeta=xsurf*dNdeta;
            wgt=WGT1(igaus);
            
            for inode=1:size(Nb,2)
                fp=mag*Nb(inode)*cross(dxdxi,dxdeta)*wgt;
                elempres([inode*3-2;inode*3-1;inode*3])=elempres([inode*3-2;inode*3-1;inode*3])+fp;
            end
        end
        
        Fpres(i0j0_surf)=Fpres(i0j0_surf)+elempres;
    end
end

