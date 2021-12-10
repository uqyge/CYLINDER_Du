
function postdisposal(nodg,Ug0,ele,nelem,nnode,ss,lambdas,DI,strain_matrix,incr)

global ELEMPRINSTRAIN
global ELEMPRINSTRESS
global ELEMSTRESS
global ELEMSTRAIN

nodg_delta=zeros(nnode,3);

ELEMPRINSTRESS=zeros(2,nelem);
ELEMPRINSTRAIN=zeros(2,nelem);
ELEMSTRESS=zeros(3,nelem);
ELEMSTRAIN=zeros(3,nelem);

h=figure(2);

for ii=1:nnode
    jj3=3*ii;
    jj2=3*ii-1;
    jj1=3*ii-2;
    nodg_delta(ii,3)=Ug0(jj3);
    nodg_delta(ii,2)=Ug0(jj2);
    nodg_delta(ii,1)=Ug0(jj1);
end

nodg_new = nodg + nodg_delta;

for ielem=1:nelem
    
    ni=ele(ielem,1);
    nj=ele(ielem,2);
    nk=ele(ielem,3);
    
    
    
    [ VV, DD ] = eig(strain_matrix(:,:,ielem));
    
    l1=VV(1,1);m1=VV(2,1);
    l2=VV(1,2);m2=VV(2,2); %主方向余弦
    
    Q = [l1^2,m1^2,l1*m1;l2^2,m2^2,l2*m2]; %坐标转换矩阵
    
    ELEMPRINSTRAIN(:,ielem) = [DD(1,1);DD(2,2)]; %主应变
    
    strain1 = max(ELEMPRINSTRAIN(:,ielem));
    strain2 = min(ELEMPRINSTRAIN(:,ielem));
    
    lambdae = [lambdas(2*ielem-1);lambdas(2*ielem)];
    
    ELEMPRINSTRESS(:,ielem) = DI*(ELEMPRINSTRAIN(:,ielem) - ss*lambdae); %主应力
    
    stress1 = max(ELEMPRINSTRESS(:,ielem));
    stress2 = min(ELEMPRINSTRESS(:,ielem));
    
    
    ELEMSTRESS(:,ielem) = Q'*ELEMPRINSTRESS(:,ielem); %应力
    ELEMSTRAIN(:,ielem) = [strain_matrix(1,1,ielem);strain_matrix(2,2,ielem);strain_matrix(1,2,ielem)];
    
    
    
    %%%%%%%%% Plot Undeformed and Deformed configurations
    
    %     plot3([nodg(ni,1),nodg(nj,1),nodg(nk,1),nodg(ni,1)],[nodg(ni,2),nodg(nj,2),nodg(nk,2),nodg(ni,2)],...
    %         [nodg(ni,3),nodg(nj,3),nodg(nk,3),nodg(ni,3)],'k-','linewidth',0.5);
    %     hold on
    
    
    plot3([nodg_new(ni,1),nodg_new(nj,1),nodg_new(nk,1),nodg_new(ni,1)],[nodg_new(ni,2),nodg_new(nj,2),nodg_new(nk,2),nodg_new(ni,2)],...
        [nodg_new(ni,3),nodg_new(nj,3),nodg_new(nk,3),nodg_new(ni,3)],'b-','linewidth',0.5);
    
    hold on
    
    if(stress1>0 && (stress2)<1e-5)
        fill3([nodg_new(ni,1),nodg_new(nj,1),nodg_new(nk,1),nodg_new(ni,1)],[nodg_new(ni,2),nodg_new(nj,2),nodg_new(nk,2),nodg_new(ni,2)],...
            [nodg_new(ni,3),nodg_new(nj,3),nodg_new(nk,3),nodg_new(ni,3)],'k');
    end
    
end


axis equal

axis off

fn = ['myfig' num2str(incr+40) '.fig'];
saveas(h,fn);

close(figure(2));

end






