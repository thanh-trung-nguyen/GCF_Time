function [Coef,x,t,Ui,Ui_tt,NoiseLevel,X_measure,freq,Use_of_Gaussian_coef] = waveeq1d_loadinput(paramfile)
% load the input parameters for the 1D wave equation

% -----------------load the parameters: 
fid = fopen(paramfile);
if fid == -1
    error('The parameter file not found');
end

text = fscanf(fid,'%s',1);  %#ok<*NASGU>
gridfile = fscanf(fid,'%s',1);

text = fscanf(fid,'%s',1); 
Xmin = str2double(fscanf(fid,'%s',1));
Xmax = str2double(fscanf(fid,'%s',1)) ;
Nx = round(str2double(fscanf(fid,'%s',1)));

text = fscanf(fid,'%s',1);
Tmax = str2double(fscanf(fid,'%s',1));
Nt = round(str2double(fscanf(fid,'%s',1)));

text = fscanf(fid,'%s',1); 
UseIncWaveForm = round(str2double(fscanf(fid,'%s',1)));

text = fscanf(fid,'%s',1); 
IncWaveFile = fscanf(fid,'%s',1);

text = fscanf(fid,'%s',1); 
freq = str2double(fscanf(fid,'%s',1));
TimeDelay = str2double(fscanf(fid,'%s',1));


text = fscanf(fid,'%s',1); 
NoiseLevel = str2double(fscanf(fid,'%s',1));

text = fscanf(fid,'%s',1); 
SolutionFile = fscanf(fid,'%s',1); 

text = fscanf(fid,'%s',1); 
X_measure = str2double(fscanf(fid,'%s',1));
SideOfExcitation = fscanf(fid,'%s',1); 

text = fscanf(fid,'%s',1); 
NrObjects = round(str2double(fscanf(fid,'%s',1)));

Xmin_obj = zeros(NrObjects,1); 
Xmax_obj = Xmin_obj;
Coef_obj = Xmin_obj;
obj_types = cell(NrObjects,1); % types of objects


for i = 1:NrObjects
    text = fscanf(fid,'%s',1); 
    Xmin_obj(i) = str2double(fscanf(fid,'%s',1));
    Xmax_obj(i) = str2double(fscanf(fid,'%s',1));
    Coef_obj(i) = str2double(fscanf(fid,'%s',1));
    obj_types{i} = fscanf(fid,'%s',1); 
end

fclose(fid);

dt = Tmax/(Nt-1); 
t = linspace(0,Tmax,Nt)';

x = linspace(Xmin,Xmax,Nx);

% test the non-uniform grid: 
% x = [linspace(0,0.25,101),linspace(0.250+0.000125,0.5,200), linspace(0.5+0.0025,1,200)];Nx = length(x);


% create the coefficient: 
if nargin < 2    
    if NrObjects >=2 
        Coef_obj = Coef_obj - 1;
        X = [Xmin_obj(1); Xmax_obj(1)];
        CoefValue = [Coef_obj(1); Coef_obj(1)];
        for i = 2:NrObjects
            if strcmpi(obj_types{i},'gaussian')
                error('The type of object must be "stepwise" for more than 1 objects');
            end
            Idx1 = find(x > Xmax_obj(i-1),1,'first');
            Idx2 = find(x < Xmin_obj(i),1,'last');
            X = [X; x(Idx1); x(Idx2); Xmin_obj(i); Xmax_obj(i)];
            CoefValue = [CoefValue; 0; 0; Coef_obj(i); Coef_obj(i)];
        end
        phi = FEM_basis(X,x);
        Coef = FEM_basis_expansion(phi,CoefValue);    
        Use_of_Gaussian_coef = 0;
    elseif NrObjects==1
        if strcmpi(obj_types{1},'gaussian')
            Coef = Gaussian_coefficient(x,Xmin_obj(1),Coef_obj(1),Xmax_obj(1));
            Use_of_Gaussian_coef = 1;
        elseif strcmpi(obj_types{1},'sine')
            
            
            
        else
            Coef_obj = Coef_obj - 1;
            X = [Xmin_obj(1); Xmax_obj(1)];
            CoefValue = [Coef_obj(1); Coef_obj(1)];
            phi = FEM_basis(X,x);
            Coef = FEM_basis_expansion(phi,CoefValue);    
            Use_of_Gaussian_coef = 0;
        end
    else % no object:       
        Coef = ones(Nx,1);
    end
    
   
else
    Coef = dlmread(CoefFile);
end

% the incident wave:
if (~UseIncWaveForm) &&  exist(IncWaveFile,'file') % load the data for incident wave if exists
    [Ui,Ui_t,Ui_tt] = load_incident_wave_from_file(IncWaveFile);
else
    [Ui,Ui_t, Ui_tt] = incident_wave(freq,x,t,SideOfExcitation,TimeDelay);
    write_incident_wave_to_file(IncWaveFile,Ui,Ui_t, Ui_tt);
end
