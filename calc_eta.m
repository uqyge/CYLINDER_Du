function [ eta ] = calc_eta( theta )
% eta function
    angle = norm(theta,2);
    if(angle<1E-10)
        eta=1/12+(1/720)*angle^2+(1/30240)*angle^4+(1/1209600)*angle^6;
    else
        eta=(1-0.5*angle*cot(0.5*angle))/(angle^2);
    end
    
end

