function output = taup_deghost(p_control,vwater,d_src,d_rec,recover_eps,nt,dt,t,omega,np,dp,p,input)

if (mod(np,2)==1)
  middlep_ind = floor(np/2)+1;
else
  middlep_ind = floor(np/2)+1;
end

tmax = t(nt);
dt_check = dt/10.0;
t_check = (0:dt_check:tmax)';
nt_check = length(t_check);


if (d_src>d_rec) 
  d_min = d_rec;
  d_max = d_src;
else
  d_min = d_src;
  d_max = d_rec;
end

zmin_check = 0.8*d_min;
zmax_check = 1.5*d_max;

taumin_check = 2.0*zmin_check/vwater;
taumax_check = 2.0*zmax_check/vwater;

ntcheck_min = floor(taumin_check/dt_check)+1;
ntcheck_max = floor(taumax_check/dt_check)+1;


d_array = zeros(7,1);



finput = fft(input);
% foutput = zeros(size(finput));
foutput = finput;


pmax = p(np);


for i_iter = 1:7
    
  index = middlep_ind+i_iter-4;
  current_p = p(index);
  
  if (p_control-(current_p*vwater)^2>0) 
  
      
    %% Check tau shift
    data = input(:,index);
    fdata = finput(:,index);
    maxvalue = zeros(nt_check,1);
    testmean = zeros(nt_check,1);
    check = zeros(nt_check,1);
    
    average = mean(abs(data));
    
    for iter = ntcheck_min:ntcheck_max
      tau = t_check(iter);
      minustshift = exp(1i*tau*omega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(omega*tau);
      denominator = denominator + 2.0;
      fP =( fdata.*F2)./denominator;
      for i = nt:-1:floor(nt/2)
        fP(i) = conj(fP(nt - i + 2));
      end

      P2 = real(ifft(fP));
      G2 = data - P2;
 
       
      maxvalue(iter) = max(P2)/average;
      test = G2 + circshift(P2,[iter-1 0]);
      testmean(iter) = average/mean(abs(test));
      check(iter) = maxvalue(iter)+testmean(iter);

    end

  
    [tmp ind] = max(check);

  
    %% Source ghost
    tau = (ind -1)*dt_check;
    d_array(i_iter) = 0.5*tau*vwater/sqrt(1-(current_p*vwater)^2);
    
        
  end
  
end


d_array(4) = 0.0;
d1 = sum(d_array)/6.0
d1 = d_src;


for i_iter = 1:np
    
  current_p = p(i_iter);
  
  if (p_control-(current_p*vwater)^2>0) 
  
    data = input(:,i_iter);
    fdata = finput(:,i_iter);
    ftemp = zeros(size(fdata));

    tau = 2.0*d1*sqrt(1-(current_p*vwater)^2)/vwater;
  
    for j_iter = 1:30
    
      minustshift = exp(1i*tau*omega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(omega*tau);
      denominator = denominator + recover_eps;
      ftmp =( fdata.*F2 + recover_eps*ftemp)./denominator;
      for i = nt:-1:floor(nt/2)
        ftmp(i) = conj(ftmp(nt - i + 2));
      end
    
    end
    
    foutput(:,i_iter) = ftmp;
  
  end
  
end


output = real(ifft(foutput));


for i_iter = 1:7
    
  index = middlep_ind+i_iter-4;
  current_p = p(index);
  
  if (p_control-(current_p*vwater)^2>0) 
  
      
    %% Check tau shift
    data = output(:,index);
    fdata = foutput(:,index);
    maxvalue = zeros(nt_check,1);
    testmean = zeros(nt_check,1);
    check = zeros(nt_check,1);
    
    
    average = mean(abs(data));
    
    for iter = ntcheck_min:ntcheck_max
      tau = t_check(iter);
      minustshift = exp(1i*tau*omega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(omega*tau);
      denominator = denominator + 2.0;
      fP =( fdata.*F2)./denominator;
      for i = nt:-1:floor(nt/2)
        fP(i) = conj(fP(nt - i + 2));
      end

      P2 = real(ifft(fP));
      G2 = data - P2;
 
      maxvalue(iter) = max(P2)/average;
      test = G2 + circshift(P2,[iter-1 0]);
      testmean(iter) = average/mean(abs(test));
      check(iter) = maxvalue(iter)+testmean(iter);

    end

  
    [tmp ind] = max(check);

  
    %% Source ghost
    tau = (ind -1)*dt_check;
    d_array(i_iter) = 0.5*tau*vwater/sqrt(1-(current_p*vwater)^2);
    
        
  end
  
end


d_array(4) = 0.0;
d2 = sum(d_array)/6.0
d2 = d_rec;



for i_iter = 1:np
    
  current_p = p(i_iter);
  
  if (p_control-(current_p*vwater)^2>0) 
  
    data = output(:,i_iter);
    fdata = foutput(:,i_iter);
    ftemp = zeros(size(fdata));

    tau = 2.0*d2*sqrt(1-(current_p*vwater)^2)/vwater;
  
    for j_iter = 1:30
    
      minustshift = exp(1i*tau*omega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(omega*tau);
      denominator = denominator + recover_eps;
      ftmp =( fdata.*F2 + recover_eps*ftemp)./denominator;
      for i = nt:-1:floor(nt/2)
        ftmp(i) = conj(ftmp(nt - i + 2));
      end
    
    end
    
    foutput(:,i_iter) = ftmp;
  
  end
  
end


output = real(ifft(foutput));




