function [r] = randa(a,b,n,m)
% make random numbers between a and b with n by m dimensions  
r = a + (b-a).*rand(n,m);
end

