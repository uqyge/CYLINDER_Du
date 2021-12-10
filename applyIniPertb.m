function [ ] = applyIniPertb( filename )
    global NODE2
    global NODEU
    global NODE
    global U
    global PARA
    
    INI_DISP = load(filename);
    NODE2(1:3,:)=NODE2(1:3,:)+INI_DISP';
    NODEU = NODE2(1:3,:) - NODE(1:3,:);
    U = reshape(INI_DISP',[PARA.NNODE*3,1]);
    
end

