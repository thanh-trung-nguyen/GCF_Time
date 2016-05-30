function f = laguerre(n,x,x0)

% calculate the Laguerre function 
% f_n(x) = Laguerre function(x - x0). 


% if nargin < 3
%     x0 = 0;
% end
% x = x - x0;
% if (n == 0)
%     f = 1;
% elseif (n == 1)
%     f = 1 - x;
% else
% %     f = 1/n*((2*n-1-x)*Laguerre_function(n-1,x) - (n-1)*Laguerre_function(n-2,x));
% 
%     f1 = 1; 
%     f2 = 1 - x;
%     for k = 1:n-1;
%         f3 = 1/(k+1)*((2*k+1-x).*f2 - k*f1);
%         f1 = f2; 
%         f2 = f3; 
%     end
%     f = f3;
% end
% f = f.*exp(-x/2);
% f = f*(-1)^n; %to make all the coefficient positives


% use the explicit formula: 
f = laguerre_explicit(n,x,x0);

