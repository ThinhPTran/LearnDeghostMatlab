function output = taup_deghost(p_control,vwater,d_src,d_rec,recover_eps,nt,dt,t,omega,np,dp,p,input)

scale = abs(0.2-recover_eps);
finput = fft(input);
foutput = zeros(size(finput));


pmax = p(np);



for i_iter = 1:np
    
  current_p = p(i_iter);
  
  if ((p_control-(current_p*vwater)^2>0)) 
    
  
    if (d_rec>0) 
      %% Source ghost
      fdata = finput(:,i_iter);
      tau = 2*d_rec*sqrt(1-(current_p*vwater)^2)/vwater;
  
      minustshift = exp(1i*tau*omega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(omega*tau);
      denominator = denominator + recover_eps;
%     denominator = denominator + recover_eps + scale*abs(current_p)/pmax;
%     denominator = ones(size(F2));
      ftmp =( fdata.*F2)./denominator;
      for i = nt:-1:floor(nt/2)
        ftmp(i) = conj(ftmp(nt - i + 2));
      end
      foutput(:,i_iter) = ftmp;
    end
    
    if (d_src>0)
      %% Receiver ghost  
      fdata = foutput(:,i_iter);
      tau = 2.0*d_src*sqrt(1-(current_p*vwater)^2)/vwater;
  
  
      minustshift = exp(1i*tau*omega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(omega*tau);
      denominator = denominator + recover_eps;
%     denominator = denominator + recover_eps + scale*abs(current_p)/pmax;
%     denominator = ones(size(F2));
      ftmp =( fdata.*F2)./denominator;
      for i = nt:-1:floor(nt/2)
        ftmp(i) = conj(ftmp(nt - i + 2));
      end
      foutput(:,i_iter) = ftmp;
    end
      
  end
  
end


output = real(ifft(foutput));
% output = sqrt(output.*abs(input));




