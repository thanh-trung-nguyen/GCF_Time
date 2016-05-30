function [RegParL2,RegParH1,NrRun,NrMeshRefine,MeshRefineFactor,MaxIter] = loadparameters_ls(paramfile)
% load the input parameters for the Laguerre's methods for the CIP:

% -----------------load the parameters: 
fid = fopen(paramfile);
if fid == -1
    error('The parameter file not found');
end

% Regularization parameters: 
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
RegParL2 = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
RegParH1 = str2double(fscanf(fid,'%s',1));

% number of runs of the least squares minimization problem
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
NrRun = round(str2double(fscanf(fid,'%s',1)));

% number of mesh refinement in each run:
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
NrMeshRefine = round(str2double(fscanf(fid,'%s',1)));

% Mesh refinement factor:
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
MeshRefineFactor = round(str2double(fscanf(fid,'%s',1)));

% Maximum number of iterations:
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
MaxIter = round(str2double(fscanf(fid,'%s',1)));

fclose(fid);
