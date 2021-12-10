function [ ] = calc_quadrature( )
    
    global GAUS3
    global WGT3
    
    global GAUS1
    global WGT1
    
    GAUS3=zeros(2,3);
    GAUS3(:,1)=[1/6;1/6];
    GAUS3(:,2)=[2/3;1/6];
    GAUS3(:,3)=[1/6;2/3];
    
    WGT3=zeros(1,3);
    WGT3(1)=1/6;
    WGT3(2)=1/6;
    WGT3(3)=1/6;
    
    GAUS1=zeros(2,1);
    GAUS1(:,1)=[1/3;1/3];
    
    WGT1=zeros(1,1);
    WGT1(1)=1/2;
end

