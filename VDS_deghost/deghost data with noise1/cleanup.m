function output = cleanup(in, out, nt)


output = zeros(size(out));


for j_iter = 1:nt
  if (abs(out(j_iter))>2.0*abs(in(j_iter)))
    out(j_iter) = in(j_iter); 
  end
end

output = out; 






