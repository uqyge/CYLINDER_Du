function [ S ] = form_S( theta )
    S=[0 -theta(3) theta(2);
       theta(3) 0 -theta(1);
       -theta(2) theta(1) 0];
end

