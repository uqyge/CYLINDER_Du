function [ ] = assemble( )
    
    global KELIST
    global PARA
    global TLIST
    
    global K
    global I
    global J
    global S
    global I0J0
    global FINT
    global BC
    global U
    
    I=zeros(PARA.NELEM*18*18,1);
    J=zeros(PARA.NELEM*18*18,1);
    S=zeros(PARA.NELEM*18*18,1);
    
    for ielem=1:PARA.NELEM
        Kel = KELIST(ielem).M;
        Tel = TLIST(ielem).M;
        
        T=zeros(18,18);
        T(1:3,1:3)=Tel;T(4:6,4:6)=Tel;T(7:9,7:9)=Tel;
        T(10:12,10:12)=Tel;T(13:15,13:15)=Tel;T(16:18,16:18)=Tel;
        
        Kel = T'*Kel*T;
       
        for m=1:18
            for n=1:18
                posimn=(ielem-1)*324+(m-1)*18+n;
                I(posimn)=I0J0(m,ielem);
                J(posimn)=I0J0(n,ielem);
                S(posimn)=Kel(m,n);
            end
        end
        
    end
    
    K = sparse(I,J,S,PARA.NNODE*6,PARA.NNODE*6);
    
    % Penalty for BC
    for i=1:PARA.NBC
        posi=BC(1,i);
        K(posi,posi)=K(posi,posi)*(1+1E8);
    end
    
    for i=1:PARA.NBC
        posi=BC(1,i);value=BC(2,i);
        FINT(posi)=FINT(posi)+(U(posi)-value);
    end
end

