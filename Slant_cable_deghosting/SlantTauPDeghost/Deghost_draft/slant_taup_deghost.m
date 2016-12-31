function output = slant_taup_deghost(p_control,vwater,x0,x_noff,alpha,eps,nt,dt,t,omega,np,dp,p,input)

finput=fft(input);
foutput=zeros(size(finput));


for i_iter = 1:np
    
  current_p=p(i_iter);
  
  if (p_control-(current_p*vwater)^2>0) 
  
    fdata=finput(:,i_iter);
    
    tau=2.0*(x_noff+x0)*tan(alpha)*sqrt(1-(current_p*vwater)^2)/vwater;
    
    minustshift = exp(1i*tau*omega);
    F2 = 1 - minustshift;
    denominator = 2 -  2*cos(omega*tau) + eps;
    
    ftmp=(fdata.*F2)./denominator;
    
      
    for i = nt:-1:floor(nt/2)
      ftmp(i) = conj(ftmp(nt - i + 2));
    end
    
    foutput(:,i_iter) = ftmp;
  
  end
  
end


output = real(ifft(foutput));