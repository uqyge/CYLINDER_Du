function [  ] = calc_localcoord_current( )

    global NODE2
    global ELEM
    global PARA
    global NODE_LOCAL
    global TLIST
       
    for ielem=1:PARA.NELEM
        x1=NODE2(1:3,ELEM(1,ielem));x2=NODE2(1:3,ELEM(2,ielem));x3=NODE2(1:3,ELEM(3,ielem));
        r1=(x2-x1)./norm(x2-x1,2);
        r3=cross((x2-x1),(x3-x1))./norm(cross((x2-x1),(x3-x1)),2);
        r2=cross(r3,r1);
        
        TLIST(ielem).M = [r1 r2 r3];

        NODE_LOCAL(ielem).x1=TLIST(ielem).M'*(x1-x1);
        NODE_LOCAL(ielem).x2=TLIST(ielem).M'*(x2-x1);
        NODE_LOCAL(ielem).x3=TLIST(ielem).M'*(x3-x1);
    end
    
end

