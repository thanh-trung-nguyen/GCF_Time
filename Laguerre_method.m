function [c,qnr,fval] = Laguerre_method(qn1,qn2,x,s,Nq,lambda,nu,options,Q_0)
% solve the inverse problem using the Laguerre's functions. 
%
% solve the Laguerre's problem in [a,b].
% qn1: Dirichlet data for q at x = b, a row vector
% qn2: Dirichlet data at x = b - h, a row vector
% x: include two ends a and b. 

if size(qn1,1) > 1
    qn1 = qn1';
end
if size(qn2,1) > 1
    qn2 = qn2';
end
if size(x,2) > 1
    x = x';
end


h = x(2)-x(1);
x = x(2:end); % remove the left end since we use the Neumann b.c. here.
Nx = length(x) - 2; % Nx is the number of unknown grid points!!!

NQ = Nx*Nq; % number of unknowns in total.

s0 = s(1); % the lower bound of the pseudo-frequency s

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


% % TEST: explicit solution:
% qn_exp = zeros(Nx+2,Nq);
% qn_exp(Nx+2,:) = qn1; qn_exp(Nx+1,:) = qn2;
% for ii = Nx+1:-1:2
%     for ik = 1:Nq        
%         for im = 1:Nq
%             for in = 1:Nq
%                 qn_exp(ii-1,ik) = qn_exp(ii-1,ik) + F(ik,im,in)*(qn_exp(ii+1,im) - qn_exp(ii,im))*(qn_exp(ii+1,in)-qn_exp(ii,in));
%             end
%             qn_exp(ii-1,ik) = qn_exp(ii-1,ik) + G(ik,im)*h*(qn_exp(ii+1,im) - qn_exp(ii,im));
%         end
%     end
%     qn_exp(ii-1,:) = -(qn_exp(ii+1,:) - 2*qn_exp(ii,:) + qn_exp(ii-1,:));
% end



% scaling the coefficient qn: 
scaling = 1e-3; 
SF = ones(1,Nq)*scaling; % scaling factor
qn1 = qn1./SF;
qn2 = qn2./SF;

if nargin >= 9
    Q_0 = Q_0/scaling;
end

% % compute the Laguerre's coefficients of the boundary data: 
Psi_n_c = zeros(1,Nq); % Neumann condition at x = c < b. FOR scattered wave it is equal to zero!


% Minimize the objective function:
alpha = 0; % Tikhonov regularization parameter:

% for storing the result:
qnr = zeros(Nx+2,Nq);

% The initial guess:   
if nargin < 9
    Q_0 = reshape(((x(1:Nx)-x(Nx+1))*qn1 - (x(1:Nx) - x(Nx+2))*qn2)/h,NQ,1); % linear initial guess
%    Q_InitGuess = reshape(ones(Nx,1)*qn2,NQ,1); % constant initial guess
end

Q_InitGuess = Q_0;


Nshift = 3; % Nshift MUST be larger than 1, only used when Niter > 1
Niter = 1;
Q1 = reshape(Q_0,Nx,Nq);
lb0 = min(Q1);
ub0 = max(Q1);

for idd = 1:Niter % layer-stripping procedure

    CarWF = (exp(-lambda*(x(Nx+1) - x(1:Nx+1)).^(nu))); % Carleman weight function


    % lower bound and upper bound of q_n:
    lb = zeros(NQ,1);
    ub = zeros(NQ,1);

    for idq = 1:Nq
        lb((idq-1)*Nx+1:idq*Nx) = lb0(idq);
        ub((idq-1)*Nx+1:idq*Nx) = ub0(idq);
    end  
    
    
% % Check the gradient calculated by the adjoint method:
% dv = 0.01; v = 8:dv:10; Nv = length(v);
% Objfun = 0*v; Grad = zeros(1,Nv); 
% Idx = 1;
% for idk = 1:Nv;
%     InitGuess = [Q_InitGuess(1:Idx-1); v(idk); Q_InitGuess(Idx+1:end)];
%     [Objfun(idk),grad] = objfun_fw(InitGuess,alpha);
%     Grad(idk) = grad(Idx);    
% end
% figure; plot(v,Objfun); title('Objective function');
% figure; plot(v,Grad(1:Nv)); hold on; plot(v(2:end-1),(Objfun(3:end)-Objfun(1:end-2))/2/dv,'--r');
% legend('Adjoint method','FD approximation'); grid on;



    
    if nargin < 9
        Q_InitGuess = lb;
    end

    % Solve the optimization problem:
    options = optimset(options,'GradObj','on','display','iter','Algorithm','sqp');
%     [Sol,fval] = fmincon(@(Q)objfun_fw(Q,alpha),Q_InitGuess,[],[],[],[],lb,ub,[],options);
    [Sol,fval] = fminunc(@(Q)objfun_fw(Q,alpha),Q_InitGuess,options);

%     % TEST: using the method for nonlinear least squares:
%     options = optimset(options,'Gradient','on','display','iter');
%     [Sol,fval] = lsqnonlin(@(Q)objfun_fw_ls(Q,alpha),Q_InitGuess,lb,ub,options);

    qn = reshape(Sol(:,end),Nx,Nq);    
    qnr(Nx+2,:)= qn1; 
    qnr(Nx+1,:) = qn2;
    qnr(1:Nx,:) = qn;
       
    % for the next layer stripping step:
    qn1 = qn(end-Nshift+2,:); qn2 = qn(end-Nshift+1,:);  
%     x = x(1:Nx-Nshift+1);   
    Nx = Nx - Nshift; NQ = Nx*Nq;

    Q_InitGuess = reshape(qn(1:end-Nshift,:),NQ,1);
    
end

qnr = qnr.*(ones(Nx+Niter*Nshift+2,1)*SF);

qnr = qnr'; qnr = [qnr(:,1) qnr];
v = compute_v_from_laguerre_coef(qnr,s(1),s(1));
c = compute_coef_from_vs(v,h,s(1));

    % Vector function of the objective functional for nonlinear least
    % squares: 
    function [J,JacJ] = objfun_fw_ls(Q,alpha)
        
        
        if nargin < 2
            alpha = 0;
        end
        Scaling = 1e2; 
        JJ = zeros(Nx+1,Nq);
        QQ = reshape(Q,Nx,Nq); QQ = [QQ; qn2; qn1];
        
        for k = 1:Nq
            i = 1; 
                JJ(i,k) = (SF(k)*(QQ(i+1,k) - QQ(i,k)) + h*Psi_n_c(k) + sum_G(G,QQ,k,i) + sum_F(F,QQ,k,i))*CarWF(i); 
                
            for i = 2:Nx
                JJ(i,k) = (SF(k)*(QQ(i+1,k) - 2*QQ(i,k) + QQ(i-1,k)) + sum_G(G,QQ,k,i) + sum_F(F,QQ,k,i))*CarWF(i);
            end   
            i = Nx+1;
                JJ(i,k) = (SF(k)*(QQ(i+1,k) - 2*QQ(i,k) + QQ(i-1,k)) + sum_G(G,QQ,k,i) + sum_F(F,QQ,k,i))*CarWF(i);

        end
        J = reshape(JJ,1,(Nx+1)*Nq)*Scaling;
        
        if nargout > 1
            JacJ = zeros((Nx+1)*Nq,NQ); 
            for k = 1:Nq
                i = 1; idx = i + (k-1)*Nx;                    
                    for l = 1:Nq
                        idl = i + (l-1)*Nx;
                        JacJ(idx,idl+1) = (h*G(k,l)*SF(l) + sum_F_grad(F,QQ,k,l,i))*CarWF(i);
                        JacJ(idx,idl) = -JacJ(idx,idl+1);
                    end
                    JacJ(idx,idx+1) = JacJ(idx,idx+1) + SF(k)*CarWF(i); 
                    JacJ(idx,idx) = JacJ(idx,idx) - SF(k)*CarWF(i);

               for i = 2:Nx-1
                    idx = i + (k-1)*Nx;
                    for l = 1:Nq
                        idl = i + (l-1)*Nx;
                        JacJ(idx,idl+1) = (h*G(k,l)*SF(l) + sum_F_grad(F,QQ,k,l,i))*CarWF(i);
                        JacJ(idx,idl) = -JacJ(idx,idl+1);
                    end
                    JacJ(idx,idx+1) = JacJ(idx,idx+1) + SF(k)*CarWF(i); 
                    JacJ(idx,idx) = JacJ(idx,idx) - 2*SF(k)*CarWF(i);
                    JacJ(idx,idx-1)   =  SF(k)*CarWF(i);                    
               end
               i = Nx; 
                   idx = i + (k-1)*Nx;
                    for l = 1:Nq
                        idl = i + (l-1)*Nx;
                        JacJ(idx,idl) = -(h*G(k,l)*SF(l) + sum_F_grad(F,QQ,k,l,i))*CarWF(i);
                    end
                    JacJ(idx,idx) = JacJ(idx,idx) - 2*SF(k)*CarWF(i);
                    JacJ(idx,idx-1)   =  SF(k)*CarWF(i);
                    
               i = Nx + 1; 
                   idx = i + (k-1)*Nx;
                    JacJ(idx,idx-1)   =  SF(k)*CarWF(i);
                    
            end
            JacJ = JacJ*Scaling;
        end
        
    end


    % objective function using the forward finite difference scheme, with
    % regularization term:
    function [J,gradJ] = objfun_fw(Q,alpha)
        
        if nargin < 2
            alpha = 0;
        end
        Scaling = 1e4; 
        J = 0; JJ = zeros(Nq,Nx+1);
        QQ = reshape(Q,Nx,Nq); QQ = [QQ; qn2; qn1];
        
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
        J = J*Scaling;
        
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
            gradJ = gradJ*Scaling;
        end
        
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


%     function stop = outfun(x,optimValues,state) % 
%         stop = false;
%  
%         switch state
%             case 'init'
%                 hold on
%             case 'iter'
%                 % Concatenate current point and objective function
%                 % value with history. x must be a row vector.
%                 history.fval = [history.fval; optimValues.fval]; % 
%                 history.x = [history.x; x'];
%             case 'done'
%                 hold off
%             otherwise
%         end
%     end

end

