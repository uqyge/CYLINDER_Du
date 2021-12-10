function [ mu ] = calc_mu( theta )

    angle=norm(theta,2);
    if(angle<1E-10)
        mu = 1/360+(1/7560)*angle^2+(1/201600)*angle^4+(1/5987520)*angle^6;
    else
        mu = (angle^2+4*cos(angle)+angle*sin(angle)-4)/(4*angle^4*sin(angle/2)^2);
    end
   
end

