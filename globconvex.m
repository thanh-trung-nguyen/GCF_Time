function [c,qnr,fval] = globconvex(Phi_n,Qn1,x,s,Nq,lambda,nu,Q0)
% solve the inverse problem using the Laguerre's functions. 


% phi_s: Dirichlet data for q at x = b. 
% Qn1: Dirichlet data at x = b - h;
% x: include c and b


h = x(2)-x(1);
x = x(2:end); 
Nx = length(x) - 2; %Nx is the number of unknown grid points!!!

NQ = Nx*Nq; % number of unknowns in total.

s0 = s(1); % the lower bound of the pseudo-frequency s
% ds2 = 0.0001; s2 = s0:ds2:1000; % for computing int_s0^infty 1/s^2*fn(s) ds

% compute Fkmn and Gkn:

filenameF = sprintf('%s%d%s%3.2f%s','globconvex_coef_F_Nq_',Nq,'_Smin_',s0,'.mat');
if ~exist(filenameF,'file')
    F = globconvex_coef_F(Nq,s0);
    eval([' save ' filenameF ' F']);
else
    eval(['load ' filenameF]);
end
filenameG = sprintf('%s%d%s%3.2f%s','globconvex_coef_G_Nq_',Nq,'_Smin_',s0,'.mat');
if ~exist(filenameG,'file')
    G = globconvex_coef_G(Nq,s0);
    eval([' save ' filenameG ' G']);
else
    eval(['load ' filenameG]);
end

% scaling the coefficient qn: 

SF = ones(1,Nq)*1e-3; % scaling factor
Phi_n = Phi_n./SF;
Qn1 = Qn1./SF;


% % compute the Laguerre's coefficients of the boundary data: 
% Phi_n = zeros(1,Nq); % Dirichlet condition at x = b
% Psi_n_b = zeros(1,Nq); % Neumann condition at x = b
Psi_n_c = zeros(1,Nq); % Neumann condition at x = c < b. FOR scattered wave it is equal to zero!

% for id = 1:Nq
% %     Phi_n(id) = laguerre_coefficient(id-1,phi_s,s(1:end-1),s0);
% %     Psi_n_b(id) = laguerre_coefficient(id-1,psi_s,s(1:end-1),s0);
%     Psi_n_c(id) = sum(laguerre(id-1,s2,s0)./(s2.^2))*ds2;
% end

% Minimize the objective function:
alpha = 0;
% Qn1 = Phi_n - h*Psi_n_b; % value of q^{N-1}_k

% for storing the result:
qnr = zeros(Nx+2,Nq);

CarWF = (exp(-lambda*(x(Nx+1) -  x(1:Nx+1)).^(nu))).^2; % Carleman weight function

for idd = 1:1 % layer-stripping procedure

    % Test the linear initial guess:    
    Q0 = reshape(((x(1:Nx)-x(Nx+1))*Phi_n - (x(1:Nx) - x(Nx+2))*Qn1)/h,NQ,1);


    % lower bound and upper bound of q_n:
    lb = zeros(NQ,1);
    ub = zeros(NQ,1);

    for idq = 1:Nq
        lb((idq-1)*Nx+1:idq*Nx) = Phi_n(idq)*ones(Nx,1);
        ub((idq-1)*Nx+1:idq*Nx) = Q0((idq-1)*Nx+1)*ones(Nx,1)*2;
    end

% %  start from the lower bound: 
%     Q0 = lb;
    
    
    % Solve the optimization problem using Quasi-Newton method:
    history.x = []; history.fval = []; 
    options = optimset('GradObj','on','display','iter','MaxIter',200,'TolFun',1e-10);
    options = optimset(options,'outputfcn',@outfun);   
    
    objfun = @(Q)objfun_fw(Q,alpha);
    % [Sol,fval] = nonlinearcg(objfun,-1,x0,tol,MaxIter);

    % Sol = fminunc(objfun,Q0,options);

    fmincon(objfun,Q0,[],[],[],[],lb,ub,[],options);
        Sol = (history.x)'; %Sol = Sol(:,end);
        fval = history.fval; %fval = fval(end);
    
    qn = reshape(Sol(:,end),Nx,Nq);    
    qnr(Nx+2,:)= Phi_n; 
    qnr(Nx+1,:) = Qn1;
    qnr(1:Nx,:) = qn;
       
    % for the next layer stripping step:
    Phi_n = Qn1; Qn1 = qn(end,:);   x = x(1:Nx+1);   Nx = Nx - 1; NQ = Nx*Nq;

end

qnr = qnr.*(ones(Nx+3,1)*SF);

v = compute_v_from_laguerre_coef(qnr',s(1),s(1));
c = compute_coef_from_vs(v,h,s(1));


    % objective function using the forward finite difference scheme, with
    % regularization term:
    function [J,gradJ] = objfun_fw(Q,alpha)
        
        if nargin < 2
            alpha = 0;
        end
        J = 0; JJ = zeros(Nq,Nx+1);
        QQ = reshape(Q,Nx,Nq); QQ = [QQ; Qn1; Phi_n];
        
        for k = 1:Nq
            i = 1; 
                JJ(k,i) = SF(k)*(QQ(i+1,k) - QQ(i,k)) + h*Psi_n_c(k) + sum_G(G,QQ,k,i) + sum_F(F,QQ,k,i);    
                J = J + JJ(k,i)^2*CarWF(i) + alpha*(QQ(i+1,k) - QQ(i,k))^2; 
            for i = 2:Nx
                JJ(k,i) = SF(k)*(QQ(i+1,k) - 2*QQ(i,k) + QQ(i-1,k)) + sum_G(G,QQ,k,i) + sum_F(F,QQ,k,i);
                J = J + JJ(k,i)^2*CarWF(i) + alpha*(QQ(i+1,k) - QQ(i,k))^2; 
            end   
            i = Nx+1;
                JJ(k,i) = SF(k)*(QQ(i+1,k) - 2*QQ(i,k) + QQ(i-1,k)) + sum_G(G,QQ,k,i) + sum_F(F,QQ,k,i);
                J = J + JJ(k,i)^2*CarWF(i); 

        end
        
        if nargout > 1
            gradJ = zeros(NQ,1); 
            for k = 1:Nq
                i = 1; idx = i + (k-1)*Nx;
                    grad = 0*gradJ;
                    for l = 1:Nq
                        idl = i + (l-1)*Nx;
                        grad(idl+1) = h*G(k,l)*SF(l) + sum_F_grad(F,QQ,k,l,i);
                        grad(idl) = -grad(idl+1);
                    end
                    grad(idx+1) = grad(idx+1) + SF(k) + alpha*2*(QQ(i+1,k) - QQ(i,k)); 
                    grad(idx) = grad(idx) - SF(k) - alpha*2*(QQ(i+1,k) - QQ(i,k));

                    gradJ = gradJ + 2*JJ(k,i)*CarWF(i)*grad;
               for i = 2:Nx-1
                    idx = i + (k-1)*Nx;
                    grad = 0*gradJ;
                    for l = 1:Nq
                        idl = i + (l-1)*Nx;
                        grad(idl+1) = h*G(k,l)*SF(l) + sum_F_grad(F,QQ,k,l,i);
                        grad(idl) = -grad(idl+1);
                    end
                    grad(idx+1) = grad(idx+1) + SF(k) + alpha*2*(QQ(i+1,k) - QQ(i,k)); 
                    grad(idx) = grad(idx) - 2*SF(k) - alpha*2*(QQ(i+1,k) - QQ(i,k));
                    grad(idx-1)   =  SF(k);
                    
                    gradJ = gradJ + 2*JJ(k,i)*CarWF(i)*grad;
               end
               i = Nx; 
                   idx = i + (k-1)*Nx;
                    grad = 0*gradJ;
                    for l = 1:Nq
                        idl = i + (l-1)*Nx;
                        grad(idl) = -( h*G(k,l)*SF(l) + sum_F_grad(F,QQ,k,l,i));
                    end
                    grad(idx) = grad(idx) - 2*SF(k) - alpha*2*(QQ(i+1,k) - QQ(i,k));
                    grad(idx-1)   =  SF(k);
                    
                    gradJ = gradJ + 2*JJ(k,i)*CarWF(i)*grad;
               i = Nx + 1; 
                   idx = i + (k-1)*Nx;
                    grad = 0*gradJ;
                    grad(idx-1)   =  SF(k);
                    
                    gradJ = gradJ + 2*JJ(k,i)*CarWF(i)*grad;            
            end
        end
        Scaling = 1e4; 
        J = J*Scaling;
        gradJ = gradJ*Scaling;
        
    end

    function Su = sum_G(G,Q2,k,i)
        Su = 0;
        for n = 1:Nq
            Su = Su + h*G(k,n)*SF(n)*(Q2(i+1,n) - Q2(i,n));
        end
    end

    function Su = sum_F(F,Q2,k,i)
        Su = 0;
        for m = 1:Nq
            for n = 1:Nq
                Su = Su + F(k,m,n)*SF(n)*SF(m)*(Q2(i+1,m) - Q2(i,m))*(Q2(i+1,n) - Q2(i,n));
            end
        end
        
    end

    function Su = sum_F_grad(F,Q2,k,l,i)
        Su = 0;
        for m = 1:Nq
            Su = Su + (F(k,l,m) + F(k,m,l))*SF(l)*SF(m)*(Q2(i+1,m) - Q2(i,m));            
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

