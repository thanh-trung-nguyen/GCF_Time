% % compute the solution of the forward problem:
% u = waveeq1d('parameter_forprob.dat','Dirichlet','coefficient_exact.dat');
[u,ce] = waveeq1d('parameter_forprob.dat','Dirichlet');

dt = 1e-4;
x =( 0:0.002:1)'; Nx= length(x); dx = x(2)-x(1);


s = 3:0.01:20; Ns = length(s);
ds = s(2) - s(1);

% compute v and q:
omega = 20;
vi = zeros(Ns,Nx);
qi = zeros(Ns,Nx);

EXP = exp(-s*2*pi/omega);
v1 = 1 - EXP;
v2 = omega^2 + s.^2;
Log = log(omega*v1./v2);
s2  = s.^2; s3 = s.^3;

v3 = 2*Log./s3 + v2./(s2.*v1).*(EXP*2*pi./omega./v2 - 2*s.*v1./(v2.^2));

for idx = 1:Nx
    vi(:,idx) = (x(idx) - x(Nx))./s + Log./s2;
    qi(:,idx) = -(x(idx) - x(Nx))./s2 - v3;  % not used!
    
end

v = compute_v(u,dt,s);
vs = v - vi; % function vs corresponds to the scattered wave

q = (vs(2:end,:) - vs(1:end-1,:))/ds;

b = 0.72; % measurement point
a = 0.5; % start of the interval in which the coefficient is estimated
nr = 2; % rate of reduce of grid points
nb = round(a/(dx*nr)) + 1; % starting index

Nq =9;
qn = laguerre_coefficient(0:Nq-1,q,s(1:end-1),s(1)); % "exact" Laguerre's coefficients




% Test the inversion: 
x2 = x(x-b < 1e-10);

x4 = x2(1:nr:end); % reduce the number of grid points
x4 = x4(nb:end); % choose a smaller domain for integration

dx4 = x4(2) - x4(1);

Nx21 = find(abs(x-x4(end)) < 1e-10);
Nx22 = find(abs(x - x4(end-1)) < 1e-10);

phi_s = q(:,Nx21);
psi_s = (q(:,Nx21) - q(:,Nx22))/dx4;

lambda = 2; % coefficient in the Carleman weight function 
nu = 1; 

%Nx2 = length(x4)- 3; % exclude one left end and two right end points
Nx2 = length(x4)- 2; % exclude two right points, but keep the left end point for a subdomain integration

NQ = Nx2*Nq; % number of unknowns in total.
% Q0 = zeros(NQ,1) + 0.03;

Q0 = reshape(qn(:,2:Nx2+1)',Nx2*Nq,1);
qn = qn(:,1:nr:Nx21)'; qn = qn(nb+1:end,:);
%Q0 = reshape(qn(2:end-2,:),NQ,1);

[cr,qnr,fval] = globconvex(qn(end,:),qn(end-1,:),x4,s,Nq,lambda,nu,Q0);

cr(cr < 1)  = 1;
vr = compute_v_from_laguerre_coef(qnr',s(1),s(1));

v2 = compute_v_from_laguerre_coef(qn',s(1),s(1));
c = compute_coef_from_vs(v2,x4(2)-x4(1),s(1));
plot([c', cr']);

cre = compute_coef_from_vs(vs(1,:),dx,s(1));
figure(1); set(gca,'fontsize',15); hold off;
plot(x,ce,'linewidth',2);
hold on; plot(x(2:end-2),cre(1:end-1),'-.k','linewidth',2); 
 plot(x4(2:end-2),cr,'--r','linewidth',2); hold off
legend('Exact coefficient','Using Laguerre expansion of the exact function v','With reconstructed function v');
grid on;
xlabel('x'); ylabel('Dielectric constant \epsilon(x)');

