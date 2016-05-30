function [X,pseudofreq,Nq,lambda,nu,CoefTruncThreshold,CoefLowerBound,CoefUpperBound,FolderName,MaxIter] = Laguerre_method_loadparameters(paramfile)
% load the input parameters for the Laguerre's methods for the CIP:

% -----------------load the parameters: 
fid = fopen(paramfile);
if fid == -1
    error('The parameter file not found');
end

% the interval in which we estimate the coefficient
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
Xmin = str2double(fscanf(fid,'%s',1));
Xmax = str2double(fscanf(fid,'%s',1)) ;
Nx = round(str2double(fscanf(fid,'%s',1)));

% interval of pseudo-frequencies:
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
Smin = str2double(fscanf(fid,'%s',1));
Smax = str2double(fscanf(fid,'%s',1)) ;
Ns = round(str2double(fscanf(fid,'%s',1)));

% number of laguerre's functions: 
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
Nq = round(str2double(fscanf(fid,'%s',1)));

% Coefficients of the Carleman weight function:
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
lambda = str2double(fscanf(fid,'%s',1));
nu = str2double(fscanf(fid,'%s',1));

if nargout > 5
    % Coefficient truncation threshold:
    text = fscanf(fid,'%s',1);  %#ok<*NASGU>
    CoefTruncThreshold = str2double(fscanf(fid,'%s',1));

    % lowerbound and upper bound
    text = fscanf(fid,'%s',1);  %#ok<*NASGU>
    CoefLowerBound = str2double(fscanf(fid,'%s',1));
    CoefUpperBound = str2double(fscanf(fid,'%s',1));
    
    % folder name: 
    text = fscanf(fid,'%s',1); FolderName = fscanf(fid,'%s',1);
    
    % maximum number of iterations
    text = fscanf(fid,'%s',1); MaxIter = round(str2double(fscanf(fid,'%s',1)));
    
end


fclose(fid);

X = linspace(Xmin,Xmax,Nx);
pseudofreq = (linspace(Smin,Smax,Ns))';
