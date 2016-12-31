function output = taup_slant_dg(recover_eps,alpha,vwater,nt,dt,npprime,dpprime,in) 


pdelta=2.0*sin(alpha)


dtpprime=0.004;
nfpprime=npprime;
dfpprime=1.0/((nfpprime)*dtpprime);
fpprime=(0:dfpprime:(nfpprime-1)*dfpprime)';
omegapprime=2.0*pi*fpprime;


t_lag=pdelta*dtpprime/dpprime


nominator=1-exp(1i*t_lag*omegapprime);
denominator=2.0-2.0*cos(t_lag*omegapprime)+recover_eps;
filter=nominator./denominator;



output=zeros(size(in));


for i_iter=1:nt


%% Test auto correlation
  tmpin=in(i_iter,:)';
  
  
  ftmpin=fft(tmpin);
  ftmpout=ftmpin.*filter;
  
  
  for j_iter=2:floor(nfpprime/2)
    ftmpout(nfpprime-j_iter+2)=conj(ftmpout(j_iter));    
  end
  
  tmpout=real(ifft(ftmpout));
  
  
  output(i_iter,:)=tmpout';
  
  
end








