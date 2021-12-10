function [ Omega ] = calc_Omega( R )

      temp=R-R';
      if(norm(temp,2)<1E-16)
          Omega=zeros(3,3);
      else
          axial=[temp(3,2);temp(1,3);temp(2,1)];
          tau=0.5*norm(axial,2);
          Omega=asin(tau)/(2*tau)*temp;
      end
end

