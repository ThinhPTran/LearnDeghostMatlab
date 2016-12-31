function output = cleanup(in, out, nx, nt, dt, nfb, fmin, fmax, res)

output = zeros(size(out));

for i_iter = 1:nx

  tmp_in = in(:,i_iter);
  tmp_out = out(:,i_iter);
  
%   size(tmp_in)
%   size(tmp_out)
    
  [tspecin,tmp] = analyse_spec_fwd(tmp_in,1,nt,dt,nfb,fmin,fmax,res);
  [tspecout,tmp] = analyse_spec_fwd(tmp_out,1,nt,dt,nfb,fmin,fmax,res);

  for j_iter = 1:nt
    if (abs(sum(tspecout(j_iter,:)))>1.0*abs(sum(tspecin(j_iter,:))))
      tspecout(j_iter,:) = 0;     
    end
  end
  
  for j_iter = 1:nfb
    output(:,i_iter) = output(:,i_iter)+(tspecout(:,j_iter));
  end

end



