function Sol = waveeq1d_FDM_core(Coef,Al,Am,Au,Bl,Bm,Bu,RightHandSide)
%
% solve the discrete system, which arises from the FDM scheme for the
% 1D wave equation,
% (M+A)Un+1 = M(2Un - Un-1) - B*Un + RighHandSide
% u(0) = 0, u'(0) = 0
% Coef must be a ROW vector
% NOTE that in the input parameters, the right hand side matrix does not
% include the initial and the end values (which are not used). Each column
% is a time-dependent data at a given value of x.


if size(Coef,1) > 1
    Coef = Coef';
end

[Nt,Nx] = size(RightHandSide);
Nt = Nt + 2; % since we remove two ends in the right hand side matrix

Up = zeros(1,Nx); % solution at the previous time step
Uc = Up; % solutions at the current and next time steps
Sol = zeros(Nt,Nx);


% solve the equation using an implicit scheme: 
Am = Am + Coef;

for n = 2:Nt-1
    % the right hand side vector: 
    RHS = RightHandSide(n-1,:) + Coef.*(2*Uc - Up) + tridiagonalmatrix_multi(Bl,Bm,Bu,Uc);
    
    Un = tridiagonal_system(Al,Am,Au,RHS);
    Sol(n+1,:) = Un; Up = Uc; Uc = Un;    
end
