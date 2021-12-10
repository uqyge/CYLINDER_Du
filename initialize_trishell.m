function [s,DI,e1,e2] = initialize_trishell( )
% initialize the global memory for 3d co-rotational membrane

global ELEM
global FORCE
global MATE
global PARA
global NODE2
global NODE

global FINT
global FEXT
global U
global I
global J
global S
global CONEX
global I0J0

global NODEU
global ROTMATRIX
global ELEMPRINSTRESS
global ELEMPRINSTRAIN
global ELEMSTRESS
global ELEMSTRAIN
global FEXTL

global FILE

ELEMPRINSTRESS=zeros(2,PARA.NELEM);
ELEMPRINSTRAIN=zeros(2,PARA.NELEM);
ELEMSTRESS=zeros(3,PARA.NELEM);
ELEMSTRAIN=zeros(3,PARA.NELEM);

FINT = zeros(PARA.NNODE*3,1);
FEXT = zeros(PARA.NNODE*3,1);
CONEX = zeros(PARA.NNODE,1);
I0J0 = zeros(9,PARA.NELEM);
NODEU = zeros(3,PARA.NNODE);
U = zeros(PARA.NNODE*3,1);
I = zeros(PARA.NELEM*9*9,1);
J = zeros(PARA.NELEM*9*9,1);
S = zeros(PARA.NELEM*9*9,1);

NODE2 = NODE;

for i=1:PARA.NFORCE
    FEXT(FORCE(1,i))=FORCE(2,i);
end
for i=1:PARA.NELEM
    CONEX(ELEM(1:3,i))=CONEX(ELEM(1:3,i))+1;
end
for i=1:PARA.NELEM
    for j=1:3
        I0J0(j*3-2:j*3,i)=ELEM(j,i)*3-2:ELEM(j,i)*3;
    end
end

E=MATE.EX;MU=MATE.PRXY;

MATE.Dm = E/(1-MU^2)*[1 MU 0;
    MU 1 0;
    0  0 (1-MU)/2];

MATE.G = E/2/(1+MU);

for i=1:PARA.NNODE
    ROTMATRIX(i).M = eye(3);
end

e1 = E; v1 = MU;
e2 = E; v2 = v1*e2/e1;

if (e1>=e2)
    s= -1;
    DI = (e1/(1-v1^2))*[1,v1;v1,1];
else
    s=1;
    DI = (e2/(1-v2^2))*[1,v2;v2,1];
end

%     FEXTL = zeros(PARA.NNODE*3,1);
%     F=7;
%     N=size(FORCE,2);
%     for i=1:PARA.NFORCE
%         FEXTL(FORCE(1,i)-1)=-F/N;
%     end


global LPRES
FEXTL = zeros(PARA.NNODE*3,1);
pres = LPRES;
%% 读取LATERAL.cdb文件获取横向压力载荷
elem_FEXTL = [];
fid = fopen(FILE.SECOND,'r');
while(~feof(fid))
    line = fgetl(fid);
    % Read EBLOCK keywords
    if(size(line,2)>6 && strcmp(line(1:6),'EBLOCK'))
        A = textscan(line,'%s','delimiter',',');
        neleml = str2double(A{1,1}{5,1});
        elem_FEXTL = zeros(neleml,15);
        line = fgetl(fid);
        for i=1:neleml
            line = fgetl(fid);
            tmp = sscanf(line,'%i%i%i%i%i%i%i%i%i%i%i%i%i%i%i');
            elem_FEXTL(i,:)=tmp;
        end
    end
end
%% 压力载荷
for ielem = 1:size(elem_FEXTL,1)
    node1 = elem_FEXTL(ielem,12);
    node2 = elem_FEXTL(ielem,13);
    node3 = elem_FEXTL(ielem,14);
    X1 = NODE(1:3,node1);
    X2 = NODE(1:3,node2);
    X3 = NODE(1:3,node3);
    AREA = norm(cross((X2-X1),(X3-X1)))*0.5;
    FP = AREA*pres/3;
    FEXTL(3*node1) = FEXTL(3*node1) -FP;
    FEXTL(3*node2) = FEXTL(3*node2) -FP;
    FEXTL(3*node3) = FEXTL(3*node3) -FP;
end







end

