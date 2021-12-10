function [  ] = calc_localcoord_initial( )

    global NODE
    global ELEM
    global PARA
    global NODE_LOCAL0
    global T0LIST
    
    for ielem=1:PARA.NELEM
        X1=NODE(1:3,ELEM(1,ielem));X2=NODE(1:3,ELEM(2,ielem));X3=NODE(1:3,ELEM(3,ielem));
        r1=(X2-X1)./norm(X2-X1,2);
        r3=cross((X2-X1),(X3-X1))./norm(cross((X2-X1),(X3-X1)),2);
        r2=cross(r3,r1);
        
        T0LIST(ielem).M = [r1 r2 r3];

        NODE_LOCAL0(ielem).X1=T0LIST(ielem).M'*(X1-X1);
        NODE_LOCAL0(ielem).X2=T0LIST(ielem).M'*(X2-X1);
        NODE_LOCAL0(ielem).X3=T0LIST(ielem).M'*(X3-X1);
    end
    
end

