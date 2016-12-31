function output = taup_deghost(p_control,vwater,d_rec,recover_eps,nt,dt,t,omega,np,dp,p,input)

finput = fft(input);
% foutput = zeros(size(finput));
foutput = finput;


for i_iter = 1:np
    
  current_p = p(i_iter);
  
  if (p_control-(current_p*vwater)^2>0) 
  
    data = input(:,i_iter);
    fdata = finput(:,i_iter);

%     tau = 2.0*d_rec*sqrt(1-(current_p*vwater)^2)/vwater;
    tau=2.0*d_rec/vwater;
  
    for j_iter = 1:1
    
      minustshift = exp(1i*tau*omega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(omega*tau) + recover_eps;
      ftmp =( fdata.*F2 )./denominator;
      for i = nt:-1:floor(nt/2)
        ftmp(i) = conj(ftmp(nt - i + 2));
      end
    
    end
    
    foutput(:,i_iter) = ftmp;
  
  end
  
end


output = real(ifft(foutput));




