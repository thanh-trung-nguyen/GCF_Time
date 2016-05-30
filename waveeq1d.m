function [TotalWave,Coef,x,X_mea] = waveeq1d(paramfile,method,Coefficient)

% solve the following problem using the FDM 
% c(x) u_tt - u_xx = ftt(x,t), 0<x<c
% u_n(0,t) = -u_t(0,t); 
% u_n(c,t) = -u_t(c,t);
% u(x,0) = u_t(x,0) = 0.

if nargin < 2
    method = 'FDM';
end
if strcmpi(method,'FDM')
    FDM = 1; FEM = 0;
elseif strcmpi(method,'FEM')
    FDM = 0; FEM = 1;
else
    error('Method is not known. Only FEM or FDM is allowed');
end

% load the input parameters:
[coef,x,t,Ui,Ui_tt,NoiseLevel,X_mea] = waveeq1d_loadinput(paramfile);
dt = t(2)-t(1);

if nargin > 2
    if length(Coefficient) == length(x)
        coef = Coefficient;
    else
        error('The length of the input coefficient is not consistent with the parameter file');
    end
end

h = x(2:end) - x(1:end-1);
   
Nx = length(x);
for nx = 1:Nx
    Ui_tt(:,nx) = (1 - coef(nx))*Ui_tt(:,nx);
end

% solve the wave equation:
if (FDM == 1) % finite difference method
    
    [Al,Am,Au,Bl,Bm,Bu] = waveeq1d_FDM_coefmat(h,dt,'ABC','ABC');     
    Sol = waveeq1d_FDM_core(coef,Al,Am,Au,Bl,Bm,Bu,dt^2*Ui_tt(2:end-1,:));
end

% % FEM NOT WORKING NOW!!!
% if (FEM == 1) % finite element method
% 
% %     % The mass matrix: 
% %     Ml = zeros(1,Nx-1); Mm = zeros(1,Nx); 
% %     for i = 1:Nx-1
% %         Ml(i) = h(i)*(c(i) + c(i+1))/12;
% %     end
% %     Mu = Ml; % upper diagonal. NOT NEEDED but provided for clarity
% %     i = 1; Mm(i) = (c(i) + c(i+1))*h(i)/6; 
% %     for i = 2:Nx-1
% %         Mm(i) = (c(i-1) + c(i))*h(i-1)/6 + (c(i) + c(i+1))*h(i)/6;
% %     end
% %     i = Nx; Mm(i) = (c(i) + c(i-1))*h(i-1)/6; 
% % 
% % 
% %     % the stiffness matrix:
% %     Bl = zeros(1,Nx-1); Bm = zeros(1,Nx); 
% %     for i = 1:Nx-1
% %         Bl(i) = -1/h(i);
% %     end
% %     Bu = Bl; % upper diagonal, not needed but provided for clarity of the code
% %     i = 1; Bm(i) = 1/h(i);
% %     for i = 2:Nx-1
% %         Bm(i) = 1/h(i-1) + 1/h(i);
% %     end
% %     i = Nx; Bm(i) = 1/h(i-1);
% % 
% %     Bl = Bl*dt^2/2; Bm = Bm*dt^2/2; Bu = Bu*dt^2/2; % for an implicit scheme
% % 
% %     % solve the equation using an implicit scheme: 
% % 
% %     Up = zeros(Nx,1); % solution at the previous time instant
% %     Uc = Up; % solution at the current time instant
% %     Sol = zeros(Nx,Nt); 
% % 
% %     Al = Ml + Bl; % lower diagonal
% %     Au = Al; % upper diagonal
% %     Am = Mm + Bm; % main diagonal
% %     Am(1) = Am(1) + dt; 
% %     Am(Nx) = Am(Nx) + dt;
% % 
% %     for n = 2:Nt-1
% %         % the right hand side vector: 
% %         RHS = zeros(Nx,1);
% %         i = 1; RHS(i) = 1/3*h(i)*ftt(i,n) + 1/6*h(i)*ftt(i+1,n);
% %         for i = 2:Nx-1
% %             RHS(i) = 1/6*h(i-1)*ftt(i-1,n) + 1/3*(h(i) + h(i-1))*ftt(i,n) + 1/6*h(i)*ftt(i+1,n);
% %         end
% %         i = Nx; RHS(i) = 1/6*h(i-1)*ftt(i-1,n) + 1/3*h(i-1)*ftt(i,n);
% % 
% %         RHS = dt^2*RHS + tridiagonalmatrix_multi(Ml,Mm,Mu,2*Uc - Up) - tridiagonalmatrix_multi(Bl,Bm,Bu,Uc);
% %         RHS(1) = RHS(1) + dt*Uc(1);
% %         RHS(Nx) = RHS(Nx) + dt*Uc(Nx);
% % 
% %         % solve the tri-diagonal system for Un: 
% %         Un = tridiagonal_system(Al,Am,Au,RHS);
% %         Sol(:,n+1) = Un; Up = Uc; Uc = Un;     
% % 
% %     end    
% end

% noiselevel = sqrt(max(h)*dt)*norm(Sol)*NoiseLevel; 
% Sol = Sol + noiselevel*2*(0.5 - rand(size(Sol))); % add random noise to the result

TotalWave = Sol + Ui; % add the incident wave to the solution
% % TotalWave = TotalWave.*(1 + 2*NoiseLevel*(0.5 - rand(size(TotalWave)))); % add random noise to the result

if nargout > 2
    Coef = coef;
end

