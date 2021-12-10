function [  ] = plotvtk_trishell( incr, scale )
% plot vtk file for display
% scale for displacement magnify ratio

    global NODE
    global ELEM
    global PARA
    global NODEU
    global ELEMPRINSTRAIN
    global ELEMPRINSTRESS
    global ELEMSTRAIN
    global ELEMSTRESS
    global CONEX
    
    fid = fopen(['L',num2str(incr),'.vtk'],'w');
    fprintf(fid,'# vtk DataFile Version 3.0\n');
    fprintf(fid,'One Hex element example\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'\n');
    
    nnode = PARA.NNODE;
    nelem = PARA.NELEM;
    
    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'%s\t%d\t%s\n','POINTS',nnode,'float');
    
    for i=1:nnode
        fprintf(fid,'%e\t%e\t%e\n',NODE(1,i)+scale*NODEU(1,i),NODE(2,i)+scale*NODEU(2,i),NODE(3,i)+scale*NODEU(3,i));
    end
    
    NCELL = nelem;
    CELL_SIZE = nelem*3+NCELL;
    
    fprintf(fid,'%s\t%d\t%d\n','CELLS',NCELL,CELL_SIZE);
    
    for i=1:nelem
        fprintf(fid,'%d\t',3);
        for j=1:3
            fprintf(fid,'%d\t',ELEM(j,i)-1);
        end
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'%s\t%d\n','CELL_TYPES',NCELL);
    for i=1:nelem
        if(2 == 2)
            fprintf(fid,'%d\n',5);
        end
    end 
    
    %Field Output
    fprintf(fid,'%s\t%d\n','POINT_DATA', nnode);
    
    fprintf(fid,'%s\t%s\t%s\n','VECTORS','displacement','float');
    for i=1:nnode
        fprintf(fid,'%e\t%e\t%e\n',NODEU(1,i),NODEU(2,i),NODEU(3,i));
    end
    
    %% Extrapolate principle vales
        nodestrain = zeros(3,PARA.NNODE);
        nodestress = zeros(3,PARA.NNODE);
        nodestatus = zeros(1,PARA.NNODE);
        nodeprinstrain = zeros(2,PARA.NNODE);
        nodeprinstress = zeros(2,PARA.NNODE);
        for ielem=1:PARA.NELEM
            nodestrain(1,ELEM(1:3,ielem))=nodestrain(1,ELEM(1:3,ielem))+[ELEMSTRAIN(1,ielem);ELEMSTRAIN(1,ielem);ELEMSTRAIN(1,ielem)]';
            nodestrain(2,ELEM(1:3,ielem))=nodestrain(2,ELEM(1:3,ielem))+[ELEMSTRAIN(2,ielem);ELEMSTRAIN(2,ielem);ELEMSTRAIN(2,ielem)]';
            nodestrain(3,ELEM(1:3,ielem))=nodestrain(3,ELEM(1:3,ielem))+[ELEMSTRAIN(3,ielem);ELEMSTRAIN(3,ielem);ELEMSTRAIN(3,ielem)]';
            nodestress(1,ELEM(1:3,ielem))=nodestress(1,ELEM(1:3,ielem))+[ELEMSTRESS(1,ielem);ELEMSTRESS(1,ielem);ELEMSTRESS(1,ielem)]';
            nodestress(2,ELEM(1:3,ielem))=nodestress(2,ELEM(1:3,ielem))+[ELEMSTRESS(2,ielem);ELEMSTRESS(2,ielem);ELEMSTRESS(2,ielem)]';
            nodestress(3,ELEM(1:3,ielem))=nodestress(3,ELEM(1:3,ielem))+[ELEMSTRESS(3,ielem);ELEMSTRESS(3,ielem);ELEMSTRESS(3,ielem)]';
        end
        nodestrain(1,:)=nodestrain(1,:)./CONEX';
        nodestrain(2,:)=nodestrain(2,:)./CONEX';
        nodestrain(3,:)=nodestrain(3,:)./CONEX';
        nodestress(1,:)=nodestress(1,:)./CONEX';
        nodestress(2,:)=nodestress(2,:)./CONEX';
        nodestress(3,:)=nodestress(3,:)./CONEX';
        
        for inode=1:PARA.NNODE
            strain = zeros(2,2);
            strain(1,1)=nodestrain(1,inode);
            strain(2,2)=nodestrain(2,inode);
            strain(1,2)=nodestrain(3,inode);
            strain(2,1)=nodestrain(3,inode);
            [~,DD]=eig(strain);
            nodeprinstrain(1,inode)=DD(1,1);
            nodeprinstrain(2,inode)=DD(2,2);
            
            stress = zeros(2,2);
            stress(1,1)=nodestress(1,inode);
            stress(2,2)=nodestress(2,inode);
            stress(1,2)=nodestress(3,inode);
            stress(2,1)=nodestress(3,inode);
            [~,DD]=eig(stress);
            nodeprinstress(1,inode)=DD(1,1);
            nodeprinstress(2,inode)=DD(2,2);
        end
        
        for i=1:PARA.NNODE
            if(nodeprinstrain(2,i)<0)
                nodestatus(i)=0;
            else
                nodestatus(i)=1;
            end
        end
        
        fprintf(fid,'%s\t%s\t%s\t%d\n','SCALARS','MembraneStatus','float',1);
        fprintf(fid,'%s\t%s\n','LOOKUP_TABLE','default');
        for inode=1:PARA.NNODE
            fprintf(fid,'%f\n',nodestatus(inode));
        end
        
        
        fprintf(fid,'%s\t%s\t%s\t%d\n','SCALARS','1stPrincipleStrain','float',1);
        fprintf(fid,'%s\t%s\n','LOOKUP_TABLE','default');
        for inode=1:PARA.NNODE
            fprintf(fid,'%f\n',nodeprinstrain(1,inode));
        end
    
        fprintf(fid,'%s\t%s\t%s\t%d\n','SCALARS','2ndPrincipleStrain','float',1);
        fprintf(fid,'%s\t%s\n','LOOKUP_TABLE','default');
        for inode=1:PARA.NNODE
            fprintf(fid,'%f\n',nodeprinstrain(2,inode));
        end
    
        fprintf(fid,'%s\t%s\t%s\t%d\n','SCALARS','1stPrincipleStress','float',1);
        fprintf(fid,'%s\t%s\n','LOOKUP_TABLE','default');
        for inode=1:PARA.NNODE
            fprintf(fid,'%f\n',nodeprinstress(1,inode));
        end
    
        fprintf(fid,'%s\t%s\t%s\t%d\n','SCALARS','2ndPrincipleStress','float',1);
        fprintf(fid,'%s\t%s\n','LOOKUP_TABLE','default');
        for inode=1:PARA.NNODE
            fprintf(fid,'%f\n',nodeprinstress(2,inode));
        end
    
    fclose(fid);
end

