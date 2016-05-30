function [propagated_data,prop_neuman_data] = data_propagation_1d(coef,x,t,DirData,NeuData)
% compute the data f(t) = u(0,t) from the given data g(t) = u(b,t), phi(t)
% = u_x(b,t) in the wave equation:
% 
% c(x) u_tt - u_xx = 0, a < x < b, 0 < t < T
% u(x,0) = u_t(x,0) = 0;
% u(b,t) = g(t)
% u_x(b,t) = phi(t)
% u(a,t) = f(t)
% 
% Assumption: u(x,T) = 0. 
% 
% Method: Use the finite difference method, and then solve from x=b to x=0.
%
% Input: 
%   coef: coefficient vector, a row vector
%   x: the spatial grid points, a row vector
%   t: vector of grid points in time, a column vector
%   DirData: a column vector of the Dirichet data at x = b
%   NeuData: a column vector of the Neumann data at x = b. 
%
% Output: the solution at x = 0. 


% compute a new spatial step size: to make sure that the explicit scheme is stable
Nt = length(t); dt = t(2) - t(1);

MaxCoef = max(coef);
dx = dt/2/sqrt(MaxCoef);

Nx = round(x(end)/dx) + 1;
xnew = linspace(x(1), x(end),Nx); 
dx = xnew(2)-xnew(1);
coef = linearinterpolation(coef,x,xnew); % interpolate the coefficient to the new grid


coef = coef*dx^2/dt^2;

u = zeros(Nt,Nx);


u(:,Nx) = DirData;
u(:,Nx-1) = DirData - dx*NeuData;

for idx = Nx-1:-1:2
    u(2:Nt-1,idx-1) = 2*u(2:Nt-1,idx) - u(2:Nt-1,idx+1) +  coef(idx)*(u(3:Nt,idx) - 2*u(2:Nt-1,idx) + u(1:Nt-2,idx));
end

propagated_data = u(:,1);

if nargout > 1
    prop_neuman_data = (u(:,2) - u(:,1))/dx;
end







