function [ui,ui_t,ui_tt] = load_incident_wave_from_file(IncWaveFile)

if strcmpi(IncWaveFile(end-3:end),'.mat')
    eval(['load ' IncWaveFile]);
    ui = Ui;
    ui_t = Ui_t;
    ui_tt = Ui_tt;
end

