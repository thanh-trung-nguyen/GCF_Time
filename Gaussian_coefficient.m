function coef = Gaussian_coefficient(x,x0,a,mu)

% compute the coefficient of the wave equation of the form: 
% c(x) = 1 + a*exp(-(x - x0)^2/mu^2);

coef = 1 + a*exp(-(x - x0).^2/mu^2);
