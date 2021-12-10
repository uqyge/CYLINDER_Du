function []=importcdb_shell(filename)
% Import ANSYS cdb file for 3d shell element

global NODE
global ELEM
global BC
global FORCE
global MATE
global PARA
global PRES

fid = fopen(filename,'r');

TMPBC = zeros();
TMPFORCE = zeros();
TEMP_PRES = zeros();

while(~feof(fid))
    line = fgetl(fid);
    
    % Read SFE
    if(size(line,2)>3 && strcmp(line(1:3),'SFE'))
        A = textscan(line,'%s','delimiter',',');
        TEMP_PRES(1,end+1)=str2double(A{1,1}{2,1});
        TEMP_PRES(2,end)=str2double(A{1,1}{3,1});
        line = fgetl(fid);
        A = textscan(line,'%s','delimiter',' ');
        TEMP_PRES(3,end)=str2double(A{1,1}{3,1});
        line = fgetl(fid);
        line = fgetl(fid);
        continue;
    end
    
    % Read NSUBST
    if(size(line,2)>6 && strcmp(line(1:6),'NSUBST'))
        A = textscan(line,'%s','delimiter',',');
        PARA.SUBSTEPS=str2double(A{1,1}{2,1});
    end
    
    % Read CNVTOL
    if(size(line,2)>6 && strcmp(line(1:6),'CNVTOL'))
        A = textscan(line,'%s','delimiter',',');
        PARA.EPS=str2double(A{1,1}{4,1});
    end
    
    % Read RBLOCK keywords
    if(size(line,2)>7 && strcmp(line(1:7),'RLBLOCK'))
        line = fgetl(fid);
        line = fgetl(fid);
        line = fgetl(fid);
        tmp = sscanf(line,'%i%i%e%e%e%e%e%e');
        PARA.h = tmp(3);
    end
    
    % Read NBLOCK keywords
    if(size(line,2)>6 && strcmp(line(1:6),'NBLOCK'))
        A = textscan(line,'%s','delimiter',',');
        PARA.NNODE = str2double(A{1,1}{5,1});
        NODE = zeros(6,PARA.NNODE);
        line = fgetl(fid);
        for i=1:PARA.NNODE
            line = fgetl(fid);
            tmp = sscanf(line,'%i%i%i%e%e%e%e%e%e%e');
            NODE(1:size(tmp,1)-3,i)=tmp(4:end);
        end
    end
    
    % Read EBLOCK keywords
    if(size(line,2)>6 && strcmp(line(1:6),'EBLOCK'))
        A = textscan(line,'%s','delimiter',',');
        PARA.NELEM = str2double(A{1,1}{5,1});
        ELEM = zeros(4,PARA.NELEM);
        line = fgetl(fid);
        for i=1:PARA.NELEM
            line = fgetl(fid);
            tmp = sscanf(line,'%i%i%i%i%i%i%i%i%i%i%i%i%i%i%i');
            ELEM(1:size(tmp,1)-11,i)=tmp(12:end);
        end
    end
    
    % Read boundary condition block
    if(size(line,2)>2 && strcmp(line(1:2),'D,'))
        tmp = textscan(line,'%s','delimiter',',');
        nodeid = str2double(tmp{1,1}{2,1});
        nodedof = tmp{1,1}{3,1};
        value = str2double(tmp{1,1}{4,1});
        if(strcmp(nodedof,'UX  '))
            TMPBC(1,end+1)=nodeid*3-2;
            TMPBC(2,end)=value;
        elseif(strcmp(nodedof,'UY  '))
            TMPBC(1,end+1)=nodeid*3-1;
            TMPBC(2,end)=value;
        elseif(strcmp(nodedof,'UZ  '))
            TMPBC(1,end+1)=nodeid*3;
            TMPBC(2,end)=value;
        end
    end
    
    % Read nodal force block
    if(size(line,2)>2 && strcmp(line(1:2),'F,'))
        tmp = textscan(line,'%s','delimiter',',');
        nodeid = str2double(tmp{1,1}{2,1});
        nodedof = tmp{1,1}{3,1};
        value = str2double(tmp{1,1}{4,1});
        if(strcmp(nodedof,'FX  '))
            TMPFORCE(1,end+1)=nodeid*3-2;
            TMPFORCE(2,end)=value;
        elseif(strcmp(nodedof,'FY  '))
            TMPFORCE(1,end+1)=nodeid*3-1;
            TMPFORCE(2,end)=value;
        elseif(strcmp(nodedof,'FZ  '))
            TMPFORCE(1,end+1)=nodeid*3;
            TMPFORCE(2,end)=value;            
        end
    end
    
    % Read material block
    if(size(line,2)>6 && strcmp(line(1:6),'MPDATA'))
        tmp = textscan(line,'%s','delimiter',',');
        name = tmp{1,1}{4,1};
        value = str2double(tmp{1,1}{7,1});
        if(strcmp(name,'EX  '))
            MATE.EX = value;
        elseif(strcmp(name,'NUXY'))
            MATE.NUXY = value;
        elseif(strcmp(name,'PRXY'))
            MATE.PRXY = value;
        end
    end
    
    % Count BC and FORCE number
    PARA.NBC = size(TMPBC,2)-1;
    PARA.NFORCE = size(TMPFORCE,2)-1;
    BC = TMPBC(:,2:end);
    FORCE = TMPFORCE(:,2:end);
    PRES = TEMP_PRES(:,2:end);
end


fclose(fid);

end



