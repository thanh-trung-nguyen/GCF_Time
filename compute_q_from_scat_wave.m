function [q,q_nu] = compute_q_from_scat_wave(scatwave,incwave,t,PseudoFreq)
% compute the function q = d/ds(ln(1 + ws/wi)) and its normal derivative; where ws and ws are the
% Laplace transforms of the scattered and incident waves, respectively
% Input: 
%   scatwave: scattered wave, given as a matrix with each column is a
%   time-dependent data at a given spatial point
%   t: time variable in the data
%   PseudoFreq: a vector of pseudo frequencies
%   incwave: incident wave, in the same format as the scattered wave.
% Output: 
%   q: the function q. Each column is the value of q at a given spatial
%   location. 
%   q_nu: normal derivative of q. 
% Nguyen Trung Thanh, UNCC 2014
% -------------------------------------------------------------------------
if size(PseudoFreq,2) > 1
    PseudoFreq = PseudoFreq'; % the pseudo frequency must be a column vector
end

ds = PseudoFreq(2) - PseudoFreq(1);
Nx = size(scatwave,2); % number of spatial points

ws = laplacetransform(scatwave,t,PseudoFreq); 
wi = laplacetransform(incwave,t,PseudoFreq);

w = 1 + ws./wi; w(w < eps) = eps;

vs = log(w)./((PseudoFreq*ones(1,Nx)).^2);

q = (vs(2:end,:) - vs(1:end-1,:))/ds;

if nargout > 1
    vs_nu = -2*ws./(ws+wi)./(PseudoFreq*ones(1,Nx));
    q_nu = (vs_nu(2:end,:) - vs_nu(1:end-1,:))/ds;
end
