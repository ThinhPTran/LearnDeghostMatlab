function y = model(withghost,omegawindow,dt,tau)

  maxP = max(withghost);
  tauint = fix(tau/dt);

  D = fft(withghost);
  average = mean(abs(withghost));

  [ntwindow tmpc] = size(withghost);

  denominator = 4.0 -  2*cos(omegawindow*tau);
  minustshift = exp(1i*tau*omegawindow);
  F = 1 - minustshift;
  
  fP =( D.*F)./denominator;
  for i = ntwindow:-1:floor(ntwindow/2)
    fP(i) = conj(fP(ntwindow - i + 2));
  end

  P = real(ifft(fP));
  G = withghost - P;
 
  
  y = G + circshift(P,[tauint 0]);




