function [f,ft,ftt] = incident_wave_form(omega,a,t)
% calculate the incident wave form: solution at one end
% f(t) = -A*(t-a)*exp(-(t-a)^2*omega^2)); where A = sqrt(2)*omega*exp(1/2).
% ft: derivative w.r.t. t
% ftt: second derivative w.r.t. t


A = sqrt(2)*omega*exp(1/2);
f = -A*(t-a).*exp(-(t-a).^2*omega^2);
ft =  A*exp(-(t-a).^2*omega^2).*(-1 + 2*omega^2*(t-a).^2);
ftt = A*exp(-(t-a).^2*omega^2).*(6 - 4*omega^2*(t-a).^2)*omega^2.*(t-a);


% case 2: 
% f(t) = sin(omega*t), t < 2*pi/omega; f(t) = 0 otherwise.
% Nt = length(t);
% f = zeros(1,Nt);
% 
% Nt2 = find(t < 2*pi/omega,1,'last');
% Nt1 = find(t > 0,1,'first');
% 
% f(Nt1:Nt2) = sin(t(Nt1:Nt2)*omega); 
% 
% if nargout > 1
%     ft = zeros(1,Nt);
%     ft(Nt1:Nt2) = omega*cos(t(Nt1:Nt2)*omega);
% end
% if nargout > 2
%     ftt = zeros(1,Nt);
%     ftt(Nt1:Nt2) = -omega^2*sin(t(Nt1:Nt2)*omega);
% end
