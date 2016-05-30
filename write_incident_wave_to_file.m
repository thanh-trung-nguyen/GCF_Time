function write_incident_wave_to_file(IncWaveFile,ui, ui_t, ui_tt)

if strcmpi(IncWaveFile(end-3:end),'.mat')
    Ui = ui;
    Ui_t = ui_t;
    Ui_tt = ui_tt;
    eval(['save ' IncWaveFile ' Ui Ui_t Ui_tt']);
end