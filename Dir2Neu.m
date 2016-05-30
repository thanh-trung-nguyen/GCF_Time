function un = Dir2Neu(ud,df,t,dist)
% compute the Neumann boundary condition given the Dirichlet condition of
% the solution to the 1D wave equation

%dist: distance from the Tx to the measurement point

Nt = length(t); dt  = t(2) - t(1);
Idx = round(dist/dt);

un = 0*ud;
un(Idx+1:Nt-1) = 2*df(1:Nt-Idx-1) - (ud(Idx+2:Nt) - ud(Idx+1:Nt-1))/dt;