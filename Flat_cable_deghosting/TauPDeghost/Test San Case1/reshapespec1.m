function output = reshapespec1(in, out, nx, nt, dt, nfb, fmin, fmax, res)

wmode=1; % Linear wmode=2; Log wmode=1

output = zeros(size(out));

sum_in=sum(in,2);
sum_out=sum(out,2);


%% Calculate gain parameters
[tspecin,freq_hz]=analyse_spec_fwd1(sum_in,wmode,nt,dt,nfb,fmin,fmax,res);
[tspecout,freq_hz]=analyse_spec_fwd1(sum_out,wmode,nt,dt,nfb,fmin,fmax,res);


ratio=zeros(1,nfb);

for i_iter=1:nfb
  maxin=max(abs(tspecin(:,i_iter)));
  maxout=max(abs(tspecout(:,i_iter)));
%   ratio(i_iter) = maxout/maxin;
  ratio(i_iter)=maxin/maxout;
end


%% Keep away from low freq part
for i_iter=1:nfb
  if (freq_hz(i_iter)<80) 
    if (ratio(i_iter) < 1)
      ratio(i_iter)=1.0; 
    end
  end
end



for i_iter = 1:nx
  
  tmp_out = out(:,i_iter);
  
  [tspecout,tmp]=analyse_spec_fwd1(tmp_out,wmode,nt,dt,nfb,fmin,fmax,res);
    
  for j_iter = 1:nfb
    if (ratio(j_iter)>0)
      tspecout(:,j_iter)=tspecout(:,j_iter)*ratio(j_iter); 
    end
  end
  
  tmp_tspecout=zeros(size(tspecout));

  for j_iter = 1:nfb
    tmp_tspecout(:,j_iter)=real(ifft(tspecout(:,j_iter)));
  end
  
  output(:,i_iter) = sum(tmp_tspecout,2);
  
end


figure();
plot(1:nfb,freq_hz);
legend('freq\_hz');



