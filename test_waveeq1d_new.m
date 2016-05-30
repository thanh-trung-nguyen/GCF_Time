% step <= 1: run the whole steps
% step = 2: run the Laguerre method and the local method
% step = 3: run only the local method
disp('Choose steps to run the test: ');
disp('<= 1: all (simulate data and then run both methods)');
disp('2: Laguerre method and Local method');
disp('3: Local method only');
step = input('step = ');

if ~exist('step','var')
    disp('Running test from beginning. Next time choose "step = 2" or "step = 3" before running the test to speed up');
    step = 1;
end

% Files of input parameters: 
parameterfile_Forward = 'parameter_forprob.dat';
parameterfile_Laguerre = 'parameter_inversion_Laguerre.dat';
parameterfile_Tikhonov = 'parameter_inversion_Tikhonov.dat';

% load the input parameters from files: 
[CoefExact,X_FDM,t,Ui,Ui_tt,NoiseLevel,X_mea,freq,Use_of_Gaussian_coef] = waveeq1d_loadinput(parameterfile_Forward); % Ui: incident wave, Ui_tt: second derivation w.r.t. t of the incident wave
[X_Lag,PseudoFreq,NrLagFunc,lambda,nu,Coef_Truncation_Threshold,coef_lowerbound,coef_upperbound,Folder,MaxIterLag] = Laguerre_method_loadparameters(parameterfile_Laguerre);
NrX_FDM = length(X_FDM); % number of spatial grid points in the FDM scheme 


% --------------------- Test parameters: to be changed
EPS = 1e-8; % small value for comparing two grid points
max_freq = freq + 5; % upper bound of the Fourier frequency of the incident wave.

% -----------------------
if Use_of_Gaussian_coef
    folder = ['Gaussian_',Folder,'_lambda',num2str(lambda)];
else
    folder = ['Stepwise_',Folder,'_lambda',num2str(lambda)];
end
if ~exist(folder,'dir')
    eval(['mkdir ' folder]);
end

% Step 1: Simulate the data: 
if step <= 1
    disp('Running test from beginning. Next time choose "step = 2" or "step = 3" before running the test to speed up');
    disp('Simulating the data: !!!');    

    % Solve the forward problem with the exact coefficient to create the data:
    TotalWave = waveeq1d(parameterfile_Forward,'FDM',CoefExact); % calculate the total wave
    
    dt = t(2)-t(1);
    dx_FDM = X_FDM(2) - X_FDM(1);
    ScatWave = TotalWave - Ui; % scattered wave.

    % ------Compute the exact Laguerre's coefficients: % for understanding the behavior of the Laguerre's coefficients only. 
    X_LagOld = X_Lag; % for visualization

    dxLag = X_Lag(2)  - X_Lag(1); % spatial step in the Laguerre's method. This may be different from the grid size of the forward solver
    ds = PseudoFreq(2) - PseudoFreq(1); % frequency step

    % Exact Laguerre's coefficients:
    q = compute_q_from_scat_wave(ScatWave,Ui,t,PseudoFreq);            
    qn = laguerre_coefficient(0:NrLagFunc-1,q,PseudoFreq(1:end-1),PseudoFreq(1)); % "exact" Laguerre's coefficients

    % ---- compute the coefficient epsilon(x) from the exact Laguerre's coefficient: 
    qnLag = zeros(NrLagFunc,length(X_Lag));
    for k = 1:NrLagFunc
        qnLag(k,:) = linearinterpolation(qn(k,:),X_FDM,X_Lag);
    end
    qLag = compute_q_from_laguerre_coef(qnLag,PseudoFreq(1:end-1),PseudoFreq(1));
    vLag = compute_v_from_laguerre_coef(qnLag,PseudoFreq(1),PseudoFreq(1));
    CoefLagExact = compute_coef_from_vs(vLag,dxLag,PseudoFreq(1));

    
    % --- display the function q and Laguerre's coefficients: 
    figure; set(gca,'fontsize',16);
    imagesc(q);     print('-depsc2',[folder,'/Exact_q.eps']);
    plot(qnLag(1:3,:)');   
    legend('q_0', 'q_1', 'q_2');
    print('-depsc2',[folder,'/Exact_Lag_Coef.eps']);


    % -----Prepare the data for inversion: 
    IdxMea = find(abs(X_FDM - X_mea) < dx_FDM/2); % the index at which the measurement is taken 
    if isempty(IdxMea)
        error('Measurement must be close or equal to a grid point');
    end
    
    noiselevel = sqrt(dx_FDM*dt)*norm(ScatWave)*NoiseLevel; 
    ScatWave = ScatWave + noiselevel*2*(0.5 - rand(size(ScatWave))); % add random noise to the result

    DirichletData = ScatWave(:,IdxMea); % measured Dirichlet data at the measurement point

    % denoising:
    DirichletData = denoising(DirichletData,t,max_freq);
   
    % truncate the data before the first peak: to avoid the effect of the
    % noise in the Laplace transform
    noiselevel_truncate = sqrt(dx_FDM*dt)*norm(ScatWave)*NoiseLevel;
    DirichletData = noise_truncate(DirichletData,noiselevel_truncate*2);   
    
    
    
end % end of Step 1.


% === Step 2: Running the Laguerre's method:
if step <= 2    
    disp('Running the Laguerre method: !!!');

    % --- Compute the boundary conditions for q:
    if abs(X_mea - X_Lag(end)) < EPS % right hand measurement
        IdxMea2 = find(abs(X_FDM - X_Lag(end-1)) <= dx_FDM/2,1,'last');
        Data1 = DirichletData; 
        %Data2 = -NeumannData*dxLag + Data1; 
    elseif abs(X_mea - X_Lag(1)) < EPS % left hand measurement:
        IdxMea2 = find(abs(X_FDM - X_Lag(2)) <= dx_FDM/2,1,'first');
        Data1 = DirichletData; 
        %Data2 = NeumannData*dxLag + Data1; 
    else    
        error('The interval in the Laguerre method must start from the location of measurement');
    end

    [q1,q1_derivative] = compute_q_from_scat_wave(Data1,Ui(:,IdxMea),t,PseudoFreq);
    q2 = q1 - q1_derivative*dxLag;
    %q2 = compute_q_from_scat_wave(Data2,Ui(:,IdxMea2),t,PseudoFreq); %this use Data2 above

    % --- Laguerre's coefficients of the boundary data: 
    qn1 = laguerre_coefficient(0:NrLagFunc-1,q1,PseudoFreq(1:end-1),PseudoFreq(1)); 
    qn2 = laguerre_coefficient(0:NrLagFunc-1,q2,PseudoFreq(1:end-1),PseudoFreq(1)); 

    % ---- Minimizing the Laguerre's functional: 
    options = optimset('MaxIter',MaxIterLag,'TolFun',1e-10);
    [CoefLag,qnr,fval] = Laguerre_method(qn1,qn2,X_Lag,PseudoFreq,NrLagFunc,lambda,nu,options);
    CoefLag(CoefLag <coef_lowerbound) = coef_lowerbound; CoefLag(CoefLag > coef_upperbound) = coef_upperbound;

    % Compute the coefficient in the whole interval, this will be used as the initial guess for the local method
    Index1 = find(X_FDM > X_Lag(2),1,'first');
    Index2 = find(X_FDM < X_Lag(end-1),1,'last');
    CoefLagFull = 0*CoefExact+1;
    CoefLagFull(Index1:Index2) = linearinterpolation((CoefLag-1)*1.8+1,X_Lag(2:end),X_FDM(Index1:Index2));  

    % --- display the function q and Laguerre's coefficients: 
    figure; set(gca,'fontsize',16);
    for n = 1:NrLagFunc
        plot(X_Lag, qnr(n,:),'linewidth',2); hold on; plot(X_Lag, qnLag(n,:),'--r','linewidth',2); hold off;
        legend('Comuted Laguerre coefficient','Exact Laguerre coefficient');
        xlabel('X'); ylabel('Laguerre coefficient');
        print('-depsc2',[folder,'/Lag_coef_',num2str(n),'.eps']);
    end

    % display the coefficient of Laguerre's method together with the true
    % one: TEST ------------------
    CoefLagTrunc = CoefLagFull; CoefLagTrunc(CoefLagFull < 0.75*max(CoefLagFull)) = 1; 
    figure; set(gca,'fontsize',16);
    idx = length(X_FDM):-1:1;
    hold off; plot(-X_FDM(idx),CoefExact(idx),'-r','linewidth',2); 
    hold on; plot(-X_FDM(idx),CoefLagFull(idx),'--b','linewidth',2); hold off;
    hold on; plot(-X_FDM(idx),CoefLagTrunc(idx),'-.k','linewidth',2); hold off; 
    legend('Exact coefficient','Result of Step 1','Result of Step 1 with 75% truncation');
    axis([-X_Lag(end) -X_FDM(1) min([min(CoefExact), min(CoefLagFull)])-0.5 max([max(CoefExact), max(CoefLagFull)]+2)]);
    xlabel('X'); ylabel('c(x)'); grid on;             
    print('-depsc2',[folder,'/Lag_method.eps']);
    

end
%  end of Laguerre's method.

% parameters for steps 3 and 4: 
[RegPar_local_L2,RegPar_local_H1,NrRun,NrIter,NrRefine,MaxIter] = loadparameters_ls(parameterfile_Tikhonov);
options = optimset('MaxIter',MaxIter,'TolFun',1e-6,'Algorithm','sqp'); % options for the fmin
        
if step <= 3    
    disp('Running the local method using Laguerre result: !!!');    
    
%     % --- The result of the Laguerre's method is used as  the initial guess
%     SubGrid = X_Lag; % 
%     InitGuess = [0; (CoefLag - 1); 0]; 

      % no mesh refinement: 
      SubGrid = X_FDM(Index1:Index2); 
      InitGuess = CoefLagFull(Index1:Index2) - 1;
   
    
    for RunIdx = 1:NrRun   
        for n = 1:NrIter
            % lower and upper bounds:     
            lb = 0*InitGuess-1 + coef_lowerbound; ub = lb + coef_upperbound;

            % run the local method:     
            [CoefLocal,CoefValues,fval,gradf] = waveeq1d_inverse_ls(Ui_tt,InitGuess,SubGrid,DirichletData,IdxMea,X_FDM,dt,lb,ub,options,RegPar_local_L2,RegPar_local_H1);

            % refine the mesh for the next iteration: uniform refinement
            SubGrid = linearinter(SubGrid,NrRefine);
            InitGuess = linearinter(CoefValues,NrRefine);

%             % Test: refine the mesh at the points of large gradient only:
%             SubGridOld = SubGrid; 
%             SubGrid = mesh_refinement_1d(SubGridOld,gradf,0.5); 
%             InitGuess = linearinterpolation(CoefValues,SubGridOld,SubGrid);            
%             plot(SubGrid,InitGuess,'*'); hold on; plot(X_FDM,CoefExact); hold off;
%             pause;

            % --- display the coefficient: 
            figure; set(gca,'fontsize',16);
            idx = length(X_FDM):-1:1;
            hold off; plot(-X_FDM(idx),CoefExact(idx),'-r','linewidth',2); 
            hold on; plot(-X_FDM(idx),CoefLocal(idx),'--b','linewidth',2); hold off;
            if RunIdx == 1
                hold on; plot(-X_FDM(idx),CoefLagFull(idx),'-.k','linewidth',2); hold off; 
                legend('Exact coefficient','Result of Step 2','Result of Step 1');
                axis([-X_Lag(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact), min(CoefLagFull)])-0.5 max([max(CoefLocal), max(CoefExact), max(CoefLagFull)]+2)]);
            else
                legend('Exact coefficient','Result of Step 3');
                axis([-X_Lag(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact)])-0.5 max([max(CoefLocal), max(CoefExact)]+2)]);
            end
            
            xlabel('X'); ylabel('c(x)'); grid on;             
            print('-depsc2',[folder,'/Run',num2str(RunIdx),'_Iter',num2str(n),'.eps']);
            
        end % finish each run        
        
         % ----after each run, we shrink the interval:
        CoefLocal(abs(CoefLocal-1) < Coef_Truncation_Threshold) = 1;  
        minvalue = min(CoefLocal); 
        maxvalue = max(CoefLocal);
        if abs(X_mea - X_Lag(end)) < EPS % right-hand side measurement                
            % find the right-most interval in which the coefficient is large 
            if minvalue >= 1 - Coef_Truncation_Threshold % strong target   
                disp('test');
                LeftIdx = find(CoefLocal == maxvalue,1,'first'); RightIdx = LeftIdx; 
            else % weak target
                LeftIdx = find(CoefLocal == minvalue,1,'first'); RightIdx = LeftIdx; 
            end                          
            while LeftIdx > 1 && abs(CoefLocal(LeftIdx) - 1) > Coef_Truncation_Threshold
                LeftIdx = LeftIdx - 1;
            end
%             LeftIdx = LeftIdx - 3                
            
            SubGrid = X_FDM(LeftIdx:IdxMea);
            InitGuess = CoefLocal(LeftIdx:IdxMea);
            NrIter = 1;
            
        elseif abs(X_mea - X_Lag(1)) < EPS % left hand measurement:
            RightIdx = find(abs(CoefLocal-1) > Coef_Truncation_Threshold,1,'first'); LeftIdx = RightIdx;
            while RightIdx < NrX_FDM && abs(CoefLocal(RightIdx)-1) > Coef_Truncation_Threshold
                RightIdx = RightIdx + 1;
            end
            %RightIdx = RightIdx + 5;       
            SubGrid = X_FDM(IdxMea:RightIdx);
            InitGuess = CoefLocal(IdxMea:RightIdx);
            NrIter = 1;
        end           
        
    end    
end % end of the local method


if step <= 4    
    disp('Running the local method alone: !!!');
    
   
%     % --- The result of the Laguerre's method is used as  the initial guess
%     SubGrid = X_Lag; % 
%     InitGuess = [0; 0*(CoefLag - 1); 0]; 
      SubGrid = X_FDM(Index1:Index2); 
      InitGuess = 0*CoefLagFull(Index1:Index2);

    ParBasis = FEM_basis(SubGrid,X_FDM); % the basis for parametrizing the coefficient
    CoefInitGuess = FEM_basis_expansion(ParBasis,InitGuess);
    
    for RunIdx = 1:NrRun   
        for n = 1:NrIter
            % lower and upper bounds:     
            lb = 0*InitGuess-1 + coef_lowerbound; ub = lb + coef_upperbound;

            % run the local method:     
            [CoefLocal,CoefValues,fval,gradf] = waveeq1d_inverse_ls(Ui_tt,InitGuess,SubGrid,DirichletData,IdxMea,X_FDM,dt,lb,ub,options,RegPar_local_L2,RegPar_local_H1);

            % refine the mesh for the next iteration: uniform refinement
            SubGrid = linearinter(SubGrid,NrRefine);
            InitGuess = linearinter(CoefValues,NrRefine);

%             % Test: refine the mesh at the points of large gradient only:
%             SubGridOld = SubGrid; 
%             SubGrid = mesh_refinement_1d(SubGridOld,gradf,0.85); 
%             InitGuess = linearinterpolation(CoefValues,SubGridOld,SubGrid);            
            
            % --- display the coefficient: 
            % --- display the coefficient: 
            figure; set(gca,'fontsize',16);
            idx = length(X_FDM):-1:1;
            hold off; plot(-X_FDM(idx),CoefExact(idx),'-r','linewidth',2); 
            hold on; plot(-X_FDM(idx),CoefLocal(idx),'--b','linewidth',2); hold off;
            if RunIdx == 1
                hold on; plot(-X_FDM(idx),CoefInitGuess(idx),'-.k','linewidth',2); hold off; 
                legend('Exact coefficient','The local method','Initial guess');
                axis([-X_LagOld(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact), min(CoefInitGuess)])-0.5 max([max(CoefLocal), max(CoefExact), max(CoefInitGuess)]+2)]);
            else
                legend('Exact coefficient','The local method');
                axis([-X_LagOld(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact)])-0.5 max([max(CoefLocal), max(CoefExact)]+2)]);
            end            
            xlabel('X'); ylabel('c(x)'); grid on;   

            print('-depsc2',[folder,'/Local_Run',num2str(RunIdx),'_Iter',num2str(n),'.eps']);
            
        end % finish each run        
        
         % ----after each run, we shrink the interval:
        CoefLocal(abs(CoefLocal-1) < Coef_Truncation_Threshold) = 1;  
        minvalue = min(CoefLocal); 
        maxvalue = max(CoefLocal);
        if abs(X_mea - X_Lag(end)) < EPS % right-hand side measurement                
            % find the right-most interval in which the coefficient is large 
            if minvalue >= 1 - Coef_Truncation_Threshold
                LeftIdx = find(CoefLocal == maxvalue,1,'first'); RightIdx = LeftIdx; 
            else 
                LeftIdx = find(CoefLocal == minvalue,1,'first'); RightIdx = LeftIdx; 
            end                          
            while LeftIdx > 1 && abs(CoefLocal(LeftIdx) - 1) > Coef_Truncation_Threshold
                LeftIdx = LeftIdx - 1;
            end
            %LeftIdx = LeftIdx - 5;                
            
            SubGrid = X_FDM(LeftIdx:IdxMea);
            InitGuess = CoefLocal(LeftIdx:IdxMea);
            NrIter = 1;
            
        elseif abs(X_mea - X_Lag(1)) < EPS % left hand measurement:
            RightIdx = find(abs(CoefLocal-1) > Coef_Truncation_Threshold,1,'first'); LeftIdx = RightIdx;
            while RightIdx < NrX_FDM && abs(CoefLocal(RightIdx)-1) > Coef_Truncation_Threshold
                RightIdx = RightIdx + 1;
            end
            %RightIdx = RightIdx + 5;       
            SubGrid = X_FDM(IdxMea:RightIdx);
            InitGuess = CoefLocal(IdxMea:RightIdx);
            NrIter = 1;
        end           
        
    end    
end % end of the local method



close all;

