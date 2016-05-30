function [Coef,CoefValues,fval] = waveeq1d_inverse_ls_DN_right(Ui_tt,InitGuess,SubGrid,DirDat,NeuDat,x,dt,lb,ub,options)
% solve the Coefficient identification problem for the 1D wave equation
% using the least squares method

% EPS = 1e-8;

if size(InitGuess,1) == 1
    InitGuess = InitGuess'; % InitGuess must be a column vector
end
Nx = length(x); % number of spatial grid points

h = x(2:end) - x(1:end-1); % mesh size
% IdxMea = find(abs(x - Xmea) < EPS); % the index at which the measurement is taken 
% NrUnknowns = length(InitGuess); 
[Al,Am,Au,Bl,Bm,Bu] = waveeq1d_FDM_coefmat(h,dt,'ABC','Neumann');     
ParBasis = FEM_basis(SubGrid,x); % the basis for parametrizing the coefficient
Ui_tt = Ui_tt(2:end-1,:)*dt^2;
NeuDat = 1/2*(NeuDat(2:end-1) + NeuDat(3:end))*dt^2/h(Nx-1);

% Check the gradient calculated by the adjoint method:
dv = 0.01; v = 2:dv:2.4; Nv = length(v);
Objfun = 0*v; Grad = zeros(1,Nv); 
for k = 1:Nv;
    [Objfun(k),grad] = objfun([v(k); InitGuess(2:end)]);
    Grad(k) = grad(1);    
end
figure; plot(v,Objfun); title('Objective function');
figure; plot(v,Grad(1:Nv)/2); hold on; plot(v(2:end-1),(Objfun(3:end)-Objfun(1:end-2))/2/dv,'--r'); hold off
legend('Adjoint method','FD approximation'); grid on;


% Solve the optimization problem using Quasi-Newton method:
history.x = []; history.fval = []; 
options = optimset(options,'outputfcn',@outfun,'GradObj','off','display','iter');


fmincon(@objfun,InitGuess,[],[],[],[],lb,ub,[],options);
Sol = (history.x)'; %Sol = Sol(:,end);
fval = history.fval; %fval = fval(end);

% [Sol,fval] = nonlinearcg(@objfun,-1,InitGuess,1e-6,6,[lb ub]);

CoefValues = Sol(:,end);
Coef = FEM_basis_expansion(ParBasis,CoefValues);


% Subfunctions: ===========================================================
    % objective functional:
    function [J,gradJ] = objfun(CoefValue)

        % compute the coefficient values at the FDM mesh nodes:
        coef = FEM_basis_expansion(ParBasis,CoefValue);
        RHS = Ui_tt;
        for nx = 1:Nx
            RHS(:,nx) = (1 - coef(nx))*Ui_tt(:,nx);
        end
        RHS(:,Nx) = RHS(:,Nx) + NeuDat;
        
        % solve the forward problem:         
        U = waveeq1d_FDM_core(coef,Al,Am,Au,Bl,Bm,Bu,RHS); % the scattered wave
        
        Residual = U(:,Nx) - DirDat;
        % compute the objective functional:
        J = sum(Residual.^2)/2*dt;
%         J = J + 1/2*ALPHA*sum((c(2:end) - c(1:end-1)).^2);
        
        
        % The gradient of the objective functional:
        if nargout > 1                        
            % solve the adjoint problem:
            RHSAdj = 0*RHS; RHSAdj(:,Nx) = Residual(3:end);
            eta = waveeq1d_FDM_core(coef,Au,Am,Al,Bu,Bm,Bl,RHSAdj(end:-1:1,:)); % the adjoint problem
            eta = eta(end:-1:3,:);
            
            % compute the gradient:
            W = -U(3:end,:) + 2*U(2:end-1,:) - U(1:end-2,:); 
            dJdC = sum((W - Ui_tt).*eta, 1)*dt;             
  
            gradJ = (dJdC*ParBasis)';
%             disp(sprintf('%s%10.9f','Norm of gradient: ', norm(gradJ)));
            
        end
    end

    function stop = outfun(x,optimValues,state) % 
        stop = false;

        switch state
            case 'init'
                hold on
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; optimValues.fval]; % 
                history.x = [history.x; x'];
            case 'done'
                hold off
            otherwise
        end
    end

end








