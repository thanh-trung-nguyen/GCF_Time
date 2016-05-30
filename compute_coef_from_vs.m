function Coef = compute_coef_from_vs(v,dx,s)


laplaceV = (v(3:end) - 2*v(2:end-1) + v(1:end-2))/dx^2;
nablaV = (v(3:end) - v(1:end-2))/2/dx;


Coef = laplaceV + s^2*nablaV.*nablaV + 2*s*nablaV + 1;

Coef = Coef';