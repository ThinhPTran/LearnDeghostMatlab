function output = improveeffect(in, out, nt)


output = zeros(size(in));



for j_iter = 1:nt
  if (abs(out(j_iter)) < 0.6*abs(in(j_iter)))
    out(j_iter) = 0.5*out(j_iter); 
  end
end

  
output=out; 



