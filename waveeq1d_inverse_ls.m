function [Coef,CoefValues,fval,gradf] = waveeq1d_inverse_ls(Ui_tt,InitGuess,SubGrid,MeaDat,IdxMea,x,dt,lb,ub,options,ALPHA,ALPHA2)
% solve the Coefficient identification problem for the 1D wave equation
% using the least squares method in the time domain
% Input: Ui_tt: second derivative of the incident wave w.r.t. time
% InitGuess: an initial guess of the coefficient at the grid points given in "SubGrid". 
% SubGrid: a set of mesh points at which the coefficient is assumed to be
% unknown. This mesh may be coarser than the FDM mesh. 
% MeaData: the measured Dirichlet data at the mesh points IdxMea.
% IdxMea: a vector of mesh points at which the data is measured
% x: FDM mesh
% dt: time step
% lb, ub: vectors of lower and upper bounds of the coefficient. These
% vectors should have the same dimension as "InitGuess".
% options: options for Matlab optimization algorithm
% ALPHA: regularization parameter in L2 norm
% ALPHA2: regularization parameter in H1 norm.
% 
% Output
% Coef: coefficient value at the FDM grid points
% CoefValues: coefficient value at the grid points "SubGrid". This is the
% same grid as the initial guess.
% fval: value of the objective function
% gradf: gradient of the objective functional at the reconstructed
% coefficient. This is used for adaptive mesh refinement. 
% 
%  Nguyen Trung Thanh UNCC 2014. 
% =========================================================================



if size(InitGuess,1) == 1
    InitGuess = InitGuess'; % InitGuess must be a column vector
end
Nx = length(x); % number of spatial grid points

h = x(2:end) - x(1:end-1); % mesh size
% IdxMea = find(abs(x - Xmea) < EPS); % the index at which the measurement is taken 
% NrUnknowns = length(InitGuess); 
[Al,Am,Au,Bl,Bm,Bu] = waveeq1d_FDM_coefmat(h,dt,'ABC','ABC');     
ParBasis = FEM_basis(SubGrid,x); % the basis for parametrizing the coefficient
NrMeaPoints = length(IdxMea);
Ui_tt = Ui_tt(2:end-1,:)*dt^2;

% % Check the gradient calculated by the adjoint method:
% dv = 0.01; v = 2.8:dv:3.2; Nv = length(v);
% Objfun = 0*v; Grad = zeros(1,Nv); 
% for k = 1:Nv;
%     [Objfun(k),grad] = objfun([v(k); InitGuess(2:end)]);
%     Grad(k) = grad(1);    
% end
% figure; plot(v,Objfun); title('Objective function');
% figure; plot(v,Grad(1:Nv)); hold on; plot(v(2:end-1),(Objfun(3:end)-Objfun(1:end-2))/2/dv,'--r');
% legend('Adjoint method','FD approximation'); grid on;

% ALPHA = 0; ALPHA2 = 1e-3; 
hmax = max(h);
c_initguess =  FEM_basis_expansion(ParBasis,InitGuess);

% Solve the optimization problem using Quasi-Newton method:
options = optimset(options,'GradObj','on','display','iter');
[Sol,fval] = fmincon(@objfun,InitGuess,[],[],[],[],lb,ub,[],options);


% % Test: adaptive mesh refinement in a loop: 
% for iter = 1:50
%     [Sol,fval] = fmincon(@objfun,InitGuess,[],[],[],[],lb,ub,[],options);
%     InitGuess = Sol(:,end);
%     [fval, gradf] = objfun(InitGuess);
%     idx = find(abs(gradf) > 0.7*grad

% [Sol,fval] = nonlinearcg(@objfun,-1,InitGuess,1e-6,6,[lb ub]);

% [J,GradJ] = objfun(Sol);

CoefValues = Sol(:,end);
Coef = FEM_basis_expansion(ParBasis,CoefValues);

if nargout > 3
    [fval,gradf] = objfun(CoefValues);
end

% Subfunctions: ===========================================================
    % objective functional:
    function [J,gradJ] = objfun(CoefValue)

        % compute the coefficient values at the FDM mesh nodes:
        c = FEM_basis_expansion(ParBasis,CoefValue);
        RHS = Ui_tt;
        for nx = 1:Nx
            RHS(:,nx) = (1 - c(nx))*Ui_tt(:,nx);
        end
        
        % solve the forward problem:         
        U = waveeq1d_FDM_core(c,Al,Am,Au,Bl,Bm,Bu,RHS);
        
        Residual = U(:,IdxMea) - MeaDat;
        % compute the objective functional:
        J = sum(sum(Residual.^2))/2*dt/NrMeaPoints;
        J = J + 1/2*ALPHA*hmax*(sum(c - c_initguess).^2) + 1/2*ALPHA2*hmax*sum((c(2:end) - c(1:end-1)).^2);
%         J = J + 1/2*ALPHA*sum((c(2:end) - c(1:end-1)).^2);
        
        
        % The gradient of the objective functional:
        if nargout > 1                        
            % solve the adjoint problem:
            RHSAdj = 0*RHS; RHSAdj(:,IdxMea) = Residual(3:end,:);
            eta = waveeq1d_FDM_core(c,Au,Am,Al,Bu,Bm,Bl,RHSAdj(end:-1:1,:)); % the adjoint problem
            eta = eta(end:-1:3,:);
            
            % compute the gradient:
            W = -U(3:end,:) + 2*U(2:end-1,:) - U(1:end-2,:); 
            dJdC = sum((W - Ui_tt).*eta, 1)*dt/NrMeaPoints;             
  
            gradJ = (dJdC*ParBasis)';
            gradJ = gradJ + (ALPHA*hmax*(c - c_initguess)'*ParBasis)' + (ALPHA2*hmax*(c(2:end) - c(1:end-1))'*(ParBasis(2:end,:) - ParBasis(1:end-1,:)))';
%             disp(sprintf('%s%10.9f','Norm of gradient: ', norm(gradJ)));
%             figure; plot(gradJ); title('gradient'); 
            
        end
    end
end








