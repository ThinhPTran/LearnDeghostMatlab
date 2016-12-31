function [w,tw] = ricker(f,dt)
%RICKER: Ricker wavelet of peak frequency f.
%
%  [w,tw] = ricker(f,dt);
%
%  IN   f : central freq. in Hz (f <<1/(2dt) )
%       dt: sampling ineterval in sec  
%
%  OUT  w:  the ricker wavelet
%       tw: axis
%
%
%  Example: 
%           [w,tw] = ricker(40,0.004); plot(tw,w);
%             
%
%  Author(s): M.D.Sacchi (sacchi@phys.ualberta.ca)
%  Copyright 1998-2003 SeismicLab
%  Revision: 1.2  Date: Dec/2002 
%  
%  Signal Analysis and Imaging Group (SAIG)
%  Department of Physics, UofA
%
%

 nw=2.2/f/dt;
 nw=2*floor(nw/2)+1;
 nc=floor(nw/2);
 w = zeros(nw,1);

 k=[1:1:nw]';

 alpha = (nc-k+1).*f*dt*pi;
 beta=alpha.^2;
 w = (1.-beta.*2).*exp(-beta);


  if nargout>1;
    tw = -(nc+1-[1:1:nw])*dt;
  end
