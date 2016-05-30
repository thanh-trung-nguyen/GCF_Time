function [ui,ui_t,ui_tt] = incident_wave(freq,x,t,SideOfExcitation,TimeDelay)
% compute the incident wave and its first and second derivatives w.r.t. t:
% ui, ui_t, ui_tt: each column is a time-dependent solution at one spatial
% point
% freq: frequency of incident wave
% x: locations at which the incident wave is evaluated
% t: time data
% SideOfExcitation: two possibility: left or right. 
% TimeDelay: a factor related to the form of the incident wave, indicating the time delay
% of the peak of the incident wave.


if nargin < 4
    SideOfExcitation = 'left'; % from the left
end

if nargin < 5
    TimeDelay = 0.2;    
end

if size(t,2) > 1
    t = t'; % t must be a COLUMN vector
end

Nx = length(x); Nt = length(t);

ui = zeros(Nt,Nx); ui_t = ui; ui_tt = ui;

if strcmpi(SideOfExcitation,'left')
    for nx = 1:Nx
        [ui(:,nx), ui_t(:,nx), ui_tt(:,nx)]  = incident_wave_form(freq,TimeDelay,t - x(nx));
    end
else
    for nx = 1:Nx
        [ui(:,nx), ui_t(:,nx), ui_tt(:,nx)]  = incident_wave_form(freq,TimeDelay,t - x(Nx) + x(nx));
    end
end    