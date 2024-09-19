function [u, gradu] = planewave(zks, pinfo, alpha)

   k1 = zks(1);
   k2 = zks(2);
   ca = cos(alpha);
   sa = sin(alpha);
   kstar = sqrt(k2^2 - (k1*ca)^2);
   c = k1*sa/kstar;
   R = (c-1)/(c+1);
   x = pinfo.r(1,:); x = x(:);
   y = pinfo.r(2,:); y = y(:);
   
   t1 = exp(1i*(k1*ca*x - k1*sa*y));
   t2 = R*exp(1i*(k1*ca*x + k1*sa*y));
   u = t1 + t2;
   [~, n] = size(pinfo.r);
   if nargout > 1
       gradu = zeros(2,n);
       gradu(1,:) = 1i*k1*ca*u.';
       gradu(2,:) = -1i*k1*sa*(t1 - t2).';
   end
end