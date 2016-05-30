function [f,df] = incidentwave(omega,t)
% Compute the incident wave form: 
% f(t) = sin(omega*t) for t < 2pi/omega, = 0 otherwise.
% and its derivative: 
% df = omega*cos(omega*t); 

Nt = length(t);
f = zeros(Nt,1); 

Idx = find(t <= 2*pi/omega);
f(Idx) = sin(omega*t(Idx));

if nargout > 1
    df = zeros(Nt,1);
    df(Idx) = omega*cos(omega*t(Idx)); 
end
