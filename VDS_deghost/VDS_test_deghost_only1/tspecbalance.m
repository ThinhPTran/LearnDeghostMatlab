function output = tspecbalance(input,nfb,nt,dt,ntwindow,wshift,thres) 


%% Prepare basic variables
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


nw = fix((nt - ntwindow)/wshift) + 1;
weight = wshift./ntwindow;
maxscale = 1/thres; 



output = zeros(size(input));



%% Calculate Norm
cwt_norm = zeros(nfb,1); 

for iter = 1:nfb
   cwt_norm(iter) = sqrt(mean(input(:,iter).*input(:,iter))); 
%     cwt_norm(iter) = max(abs(time_spectral(:,iter))); 
end



%% Normalize time_spectral
for iter = 1:nfb
   if (cwt_norm(iter)>0)
     input(:,iter) = input(:,iter)/cwt_norm(iter);
   end
end



%% Diversity 
for twidx = 1:nw
        
  twidx;
        
  fidx = (twidx-1)*wshift + 1;
  lidx = (twidx-1)*wshift + ntwindow;
  
  tmp_in = input(fidx:lidx,:); 
  meaninput = zeros(nfb,1); 
  scale = zeros(nfb,1); 
  
  
  % Calculate some parms
  for i_iter = 1:nfb
    meaninput(i_iter) = mean(abs(tmp_in(:,i_iter)));  
  end
  
  
  % Calculate scale
  target_pow = mean(meaninput); 
  
  for i_iter = 1:nfb
      
     if (abs(meaninput(i_iter))>0) 
       scale(i_iter) = target_pow/meaninput(i_iter);
     else 
       scale(i_iter) = 0; 
     end
     
     if (scale(i_iter) > maxscale) 
        scale(i_iter) = maxscale; 
%      elseif (scale(i_iter) < 1.0) 
%         scale(i_iter) = 1.0; 
     end
     
  end
  
  
  % Scale back
  for i_iter = 1:nfb
     tmp_in(:,i_iter) = tmp_in(:,i_iter)*scale(i_iter); 
  end
  

  for i_iter = 1:nfb 
    for j_iter = 1:ntwindow
      output(fidx + j_iter - 1,i_iter) = output(fidx + j_iter - 1,i_iter) + weight*tmp_in(j_iter,i_iter);
    end
  end
 
end



% for iter = 1:nfb
%   figure();
%   plot(t,input(:,iter),t,output(:,iter)); 
%   legend('time\_spectral','time\_spectral\_out'); 
% end



%% Renormalize
for iter = 1:nfb
   if (cwt_norm(iter)>0)
     output(:,iter) = output(:,iter)*cwt_norm(iter);
   end
end




