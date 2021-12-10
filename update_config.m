function [  ] = update_config()
% update geometry
    
    global NODE
    global NODEU
    global NODE2
    global U
    global PARA
    
    for i = 1:size(U,1)
        if(abs(U(i))<1E-20)
            U(i) = 0;
        end
    end
    
    for i=1:PARA.NNODE
        posi=3*(i-1);
        NODE2(1:3,i)=NODE(1:3,i)+U(posi+1:posi+3);
        NODEU(1:3,i)=U(posi+1:posi+3);
    end
    
end

