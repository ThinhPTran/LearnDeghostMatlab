function output = cleanup(in, out, nt, dt, nfb, fmin, fmax, res)


output = zeros(size(out));


tmp_in = in;
tmp_out = out;
  
    
[tspecin,tmp]=analyse_spec_fwd(tmp_in,1,nt,dt,nfb,fmin,fmax,res);
[tspecout,tmp]=analyse_spec_fwd(tmp_out,1,nt,dt,nfb,fmin,fmax,res);


for j_iter = 1:nt
  if (abs(sum(tspecout(j_iter,:)))>2.0*abs(sum(tspecin(j_iter,:))))
    tspecout(j_iter,:) = 0;     
  end
end

  
for j_iter = 1:nfb
  output = output+(tspecout(:,j_iter));
end




