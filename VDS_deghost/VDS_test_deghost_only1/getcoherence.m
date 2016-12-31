function [output P1out P2out P3out] = getcoherence(input, P1, P2, P3, nt, dt, nfb, fmin, fmax, res)


output = zeros(size(input));
P1out = zeros(size(input)); 
P2out = zeros(size(input)); 
P3out = zeros(size(input)); 


for j_iter = 1:nt
  if ( ...
       (input(j_iter)*P1(j_iter) < 0 ) || ...
       (input(j_iter)*P2(j_iter) < 0 ) || ...   
       (input(j_iter)*P3(j_iter) < 0 ) || ...
       (P1(j_iter)*P2(j_iter) < 0 ) || ...
       (P1(j_iter)*P3(j_iter) < 0 ) || ...
       (P2(j_iter)*P3(j_iter) < 0 ) ...
     )
    %tspecout(j_iter,:) = tspecin(j_iter,:);  
    P1(j_iter,:) = 0; 
    P2(j_iter,:) = 0; 
    P3(j_iter,:) = 0; 
  end
end

  
output=input; 
P1out = P1; 
P2out = P2; 
P3out = P3; 






