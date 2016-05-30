function [Ll,Lm,Lu,Rl,Rm,Ru] = waveeq1d_FDM_coefmat(h,dt,bdc_left,bdc_right)
% compute the coefficient matrix in the FDM for the wave equation with ABC or Neumann BC
% L: left hand side matrix
% R: right hand side matrix


if nargin < 3
    bdc_left = 'ABC';
    bdc_right = 'ABC';    
end

Nx = length(h)+1;

% create the coefficient matrix: 
Al = zeros(1,Nx-1); Am = zeros(1,Nx); Au = Al;

Am(1) = 1/h(1)^2; Au(1) = -Am(1);
for i = 2:Nx-1
    Al(i-1) = -2/(h(i-1) + h(i))/h(i-1); 
    Am(i) = 2/h(i-1)/h(i);
    Au(i) = -2/(h(i-1) + h(i))/h(i);
end
Am(Nx) = 1/h(Nx-1)^2; Al(Nx-1) = -Am(Nx);

Al = Al*dt^2/2; Am = Am*dt^2/2; Au = Au*dt^2/2;


Ll = Al; Lm = Am; Lu = Au;
Rl = -Al; Rm = -Am; Ru = -Au;

if strcmpi(bdc_left,'ABC')
    Lm(1) = Lm(1) + dt/h(1); 
    Rm(1) = Rm(1) + dt/h(1); 
elseif strcmpi(bdc_left,'Dirichlet')
    Lm(1) = Lm(1) + dt^2/h(1)^2/2; 
    Rm(1) = Rm(1) - dt^2/h(1)^2/2;
end
if strcmpi(bdc_right,'ABC')
    Lm(Nx) = Lm(Nx) + dt/h(Nx-1);
    Rm(Nx) = Rm(Nx) + dt/h(Nx-1);
elseif strcmpi(bdc_right,'Dirichlet')
    Lm(Nx) = Lm(Nx) + dt^2/h(Nx-1)^2/2;
    Rm(Nx) = Rm(Nx) - dt^2/h(Nx-1)^2/2;
end

