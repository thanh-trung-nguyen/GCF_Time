function test_ls

% calculate the measured data: 
Ce = 4; A = 10; dt = 0.001; t = 0:dt:1;
f = sin(2*pi*t); f = dt*f(2:end-1);

MeaDat = forwardsolver(Ce,A,dt,f);


dv = 0.01; v = 2:dv:6; Nv = length(v);
Objfun = 0*v; Grad = zeros(1,Nv); 
for k = 1:Nv;
    [Objfun(k),Grad(k)] = objfun(v(k));
end
figure; plot(v,Objfun); title('Objective function');
figure; plot(v,Grad(1:Nv)); hold on; plot(v(2:end-1),(Objfun(3:end)-Objfun(1:end-2))/2/dv,'--r');
legend('Adjoint method','FD approximation'); grid on;



% Subfunctions: ===========================================================
    % objective functional:
    function [J,gradJ] = objfun(C)
    
        % solve the forward problem:         
        U = forwardsolver(C,A,dt,f); % the scattered wave
        
        Residual = U - MeaDat;
        % compute the objective functional:
        J = sum(Residual.^2)/2*dt;
        
        % The gradient of the objective functional:
        if nargout > 1                        
            % solve the adjoint problem:
            RHSAdj = Residual;
            eta = forwardsolver(C,A,dt,RHSAdj(end:-1:3)); % the adjoint problem
            eta = eta(end:-1:3);
            
            % compute the gradient:
            W = -U(3:end) + 2*U(2:end-1) - U(1:end-2); 
            gradJ = sum(W.*eta)*dt;             
           
        end
    end

    function Sol = forwardsolver(C,A,dt,f)
        Nt = length(f) + 2; 
        Sol = zeros(Nt,1); 
        A = A*dt/2;
        for n = 2:Nt-1
            Sol(n+1) = (C*(2*Sol(n) - Sol(n-1)) - A*Sol(n) + f(n-1))/(C + A);
        end
    end

end








