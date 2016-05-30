function run_Laguerre_method_1d(parameterfile_Forward,parameterfile_Laguerre,parameterfile_Tikhonov)
% run the combination of the Laguerre's method and the local method for 1d
% inverse problem.
% 
%



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

% % Files of input parameters: 
% parameterfile_Forward = 'parameter_forprob.dat';
% parameterfile_Laguerre = 'parameter_inversion_Laguerre.dat';
% parameterfile_Tikhonov = 'parameter_inversion_Tikhonov.dat';

% Test parameters: 
EPS = 1e-8;
Coef_Truncation_Threshold = 0.2; % truncate the coefficient of the local method for the next run of the Laguerre's method
max_freq = 30; % upper bound of the Fourier frequency of the incident wave.

folder ='example2';
if ~exist(folder,'dir')
    eval(['mkdir ' folder]);
end

% Step 1: Simulate the data: 
if step <= 1
    disp('Running test from beginning. Next time choose "step = 2" or "step = 3" before running the test to speed up');
    disp('Simulating the data: !!!');    

    % Solve the forward problem with the exact coefficient to create the data:
    [CoefExact,X_FDM,t,Ui,Ui_tt,NoiseLevel,X_mea] = waveeq1d_loadinput(parameterfile_Forward); % Ui: incident wave, Ui_tt: second derivation w.r.t. t of the incident wave
    CoefExact = Gaussian_coefficient(X_FDM,0.58,1,0.04);
    TotalWave = waveeq1d(parameterfile_Forward,'FDM',CoefExact); % calculate the total wave
    
    dt = t(2)-t(1);
    dx_FDM = X_FDM(2) - X_FDM(1);
    ScatWave = TotalWave - Ui; % scattered wave.

    % -----Prepare the data for inversion: 
    IdxMea = find(abs(X_FDM - X_mea) < dx_FDM/2); % the index at which the measurement is taken 
    if isempty(IdxMea)
        error('Measurement must be close or equal to a grid point');
    end
    
    DirichletData = ScatWave(:,IdxMea); % measured Dirichlet data at the measurement point
    NeumannData = (ScatWave(:,IdxMea+1) - ScatWave(:,IdxMea))/(dx_FDM); % Neuman condition at the measurement point

   % denoising: 
   DirichletData = denoising(DirichletData,t,max_freq);
   
   noiselevel_truncate = sqrt(dx_FDM*dt)*norm(ScatWave)*NoiseLevel;
   DirichletData = noise_truncate(DirichletData,noiselevel_truncate*2);
   
    
%     % -----test the forward solver with Neumann b.c.:
%     coef = CoefExact(1:IdxMea);
%     
%     RHS = Ui_tt(2:end-1,1:IdxMea);
%         for nx = 1:IdxMea
%             RHS(:,nx) = (1 - coef(nx))*Ui_tt(2:end-1,nx)*dt^2;
%         end
%     h = X_FDM(2:IdxMea) - X_FDM(1:IdxMea-1);
%     [Al,Am,Au,Bl,Bm,Bu] = waveeq1d_FDM_coefmat(h,dt,'ABC','Neumann');     
%     RHS(:,IdxMea) = RHS(:,IdxMea) + 1/4*(NeumannData(2:end-1) + NeumannData(3:end))*dt^2/h(end); % Neumann data
% 
% %     [Al,Am,Au,Bl,Bm,Bu] = waveeq1d_FDM_coefmat(h,dt,'ABC','Dirichlet');     
% %     RHS(:,IdxMea) = RHS(:,IdxMea) + 1/2*(ScatWave(2:end-1,IdxMea+1) + ScatWave(3:end,IdxMea+1))*dt^2/h(end)^2; % Dirichlet data
% 
%     U = waveeq1d_FDM_core(coef,Al,Am,Au,Bl,Bm,Bu,RHS); % the scattered wave

    % ------Compute the exact Laguerre's coefficients:
    [X_Lag,PseudoFreq,NrLagFunc,lambda,nu] = Laguerre_method_loadparameters(parameterfile_Laguerre);
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

end % end of Step 1.


% === Step 2: Running the Laguerre's method AND the local method:
if step <= 2
    NrX_FDM = length(X_FDM);
    for RunIdx = 1:2 % run the Laguerre's method and the local method twice!
       
        disp('Running the Laguerre method: !!!');
        [X_Lag,PseudoFreq,NrLagFunc,lambda,nu] = Laguerre_method_loadparameters(parameterfile_Laguerre);
 
        % ----shrink the spatial interval for the second run of the Laguerre's method
        if RunIdx > 1 
            CoefLocal(abs(CoefLocal-1) < Coef_Truncation_Threshold) = 1;  %#ok<SAGROW>
            if abs(X_mea - X_Lag(end)) < EPS % right-hand side measurement
                LeftIdx = find(abs(CoefLocal-1) > Coef_Truncation_Threshold,1,'last'); RightIdx = LeftIdx; 
                while LeftIdx > 1 && abs(CoefLocal(LeftIdx) - 1) > Coef_Truncation_Threshold
                    LeftIdx = LeftIdx - 1;
                end
                LeftIdx = LeftIdx - 5;
                LeftIdxLag = find(X_Lag < X_FDM(LeftIdx),1,'last');
                dxLag = dxLag/2;
                NxLag = round((X_Lag(end)-X_Lag(LeftIdxLag))/dxLag) + 1;
                X_Lag  = linspace(X_Lag(LeftIdxLag), X_Lag(end),NxLag);
                
            elseif abs(X_mea - X_Lag(1)) < EPS % left hand measurement:
                RightIdx = find(abs(CoefLocal-1) > Coef_Truncation_Threshold,1,'first'); LeftIdx = RightIdx;
                while RightIdx < NrX_FDM && abs(CoefLocal(RightIdx)-1) > Coef_Truncation_Threshold
                    RightIdx = RightIdx + 1;
                end
                RightIdx = RightIdx + 5;
                RightIdxLag = find(X_Lag > X_FDM(RightIdx),1,'first');
                dxLag = dxLag/2;
                NxLag = round(-(X_Lag(1)+X_Lag(RightIdxLag))/dxLag) + 1;
                X_Lag  = linspace(X_Lag(1), X_Lag(RightIdxLag),NxLag);
                
            else    
                error('One end of the Laguerre method must be equal to the location of measurement');
            end
            
            % compute the function q for the reconstructed coefficient (initial guess for the second run of the Laguerre's method): 
            TotalWave2 = waveeq1d(parameterfile_Forward,'FDM',CoefLocal); % calculate the total wave
            ScatWave2 = TotalWave2 - Ui; % scattered wave.
            q_step2 = compute_q_from_scat_wave(ScatWave2,Ui,t,PseudoFreq);            
            qn = laguerre_coefficient(0:NrLagFunc-1,q_step2,PseudoFreq(1:end-1),PseudoFreq(1)); % "exact" Laguerre's coefficients
            qnLag = zeros(NrLagFunc,length(X_Lag)); % we use qnLag as the initial guess for the Laguerre's method. 
            for k = 1:NrLagFunc
                qnLag(k,:) = linearinterpolation(qn(k,:),X_FDM,X_Lag);
            end
        end
     

        % ----- Laguerre's method: compute the boundary conditions for q:
        if abs(X_mea - X_Lag(end)) < EPS % right hand measurement
            IdxMea2 = find(abs(X_FDM - X_Lag(end-1)) <= dx_FDM/2,1,'last');
            Data1 = DirichletData; 
            Data2 = -NeumannData*dxLag + Data1; 
        elseif abs(X_mea - X_Lag(1)) < EPS % left hand measurement:
            IdxMea2 = find(abs(X_FDM - X_Lag(2)) <= dx_FDM/2,1,'first');
            Data1 = DirichletData; 
            Data2 = NeumannData*dxLag + Data1; 
        else    
            error('The interval in the Laguerre method must start from the location of measurement');
        end
        
        [q1,q1_derivative] = compute_q_from_scat_wave(Data1,Ui(:,IdxMea),t,PseudoFreq);
%         q2 = compute_q_from_scat_wave(Data2,Ui(:,IdxMea2),t,PseudoFreq);
        q2 = q1 - q1_derivative*dxLag;
      
        
        qn1 = laguerre_coefficient(0:NrLagFunc-1,q1,PseudoFreq(1:end-1),PseudoFreq(1)); % Laguerre's coefficients of the boundary data
        qn2 = laguerre_coefficient(0:NrLagFunc-1,q2,PseudoFreq(1:end-1),PseudoFreq(1)); 

        % ---- Minimizing the Laguerre's functional: 
        options = optimset('MaxIter',1000,'TolFun',1e-10);
        if RunIdx == 1 
            [CoefLag1,qnr,fval] = Laguerre_method(qn1,qn2,X_Lag,PseudoFreq,NrLagFunc,lambda,nu,options);
            CoefLag = CoefLag1;
        else
            [CoefLag2,qnr,fval] = Laguerre_method(qn1,qn2,X_Lag,PseudoFreq,NrLagFunc,0,nu,options,reshape(qnLag(:,2:end-2)',(length(X_Lag)-3)*NrLagFunc,1));
            CoefLag = CoefLag2;
        end

        % Compute the coefficient in the whole interval, this will be used as the initial guess for the local method
        Index1 = find(X_FDM > X_Lag(2),1,'first');
        Index2 = find(X_FDM < X_Lag(end-1),1,'last');
        CoefLagFull = 0*CoefExact+1;
        CoefLagFull(Index1:Index2) = linearinterpolation(CoefLag,X_Lag(2:end),X_FDM(Index1:Index2));  


        % 2. Estimate the coefficient using the local method with adaptive parametrization
        disp('Running the local method:');

        % recursively adaptive parametrization:
        NrRefine = 2; % refinement factor
        
        if RunIdx==1
            NrIter = 3; % number of total refinements
%             NrUnknowns = length(CoefLag1);
            SubGrid = X_Lag(2:end-1);
            InitGuess = (CoefLag - 1); 
        else
            NrIter = 1; % number of total refinements
%             NrUnknowns = length(CoefLag2);
%             SubGrid = linspace(X_Lag(2), X_Lag(end-1),NrUnknowns);
%             InitGuess = linearinterpolation(CoefLag2,X_Lag(2:end-1),SubGrid) - 1; 
% %             InitGuess = max(CoefLag)*ones(NrUnknowns,1) - 1; % start from the homogeneous medium

            SubGrid = X_Lag(2:end-1);
            InitGuess = (CoefLag - 1);         
        end            
            
        ParBasis = FEM_basis(SubGrid,X_FDM); % the basis for parametrizing the coefficient
        CoefInitGuess = FEM_basis_expansion(ParBasis,InitGuess);

        options = optimset('MaxIter',150,'TolFun',1e-6,'Algorithm','sqp'); % options for the fmin

        for n = 1:NrIter
            fprintf('%s%d','Nr of refinements: ', n);
            lb = 0*InitGuess - 1; ub = lb + 30;

            [CoefLocal,CoefValues,fval,gradf] = waveeq1d_inverse_ls(Ui_tt,InitGuess,SubGrid,DirichletData,IdxMea,X_FDM,dt,lb,ub,options);

            % refine the mesh: uniform refinement
            SubGrid = linearinter(SubGrid,NrRefine);
            InitGuess = linearinter(CoefValues,NrRefine);

%             % Test: refine the mesh at the points of large gradient only:
%             SubGridOld = SubGrid; 
%             SubGrid = mesh_refinement_1d(SubGridOld,gradf,0.85); 
%             InitGuess = linearinterpolation(CoefValues,SubGridOld,SubGrid);            
   

            % --- display the coefficient: 
            figure; set(gca,'fontsize',16);
            idx = length(X_FDM):-1:1;
            hold off; plot(X_FDM,CoefLocal(idx),'--b','linewidth',2); 
            hold on; plot(X_FDM,CoefLagFull(idx),'-k','linewidth',2); hold off; 
            hold on; plot(X_FDM,CoefExact(idx),'-.r','linewidth',2); hold off; 

            legend('Reconstruction using the local method','Reconstruction using Laguerre method','Exact coefficient');
            axis([X_FDM(1) X_FDM(end) min([min(CoefLocal), min(CoefExact), min(CoefInitGuess), min(CoefLagFull)])-0.5 max([max(CoefLocal), max(CoefExact), max(CoefInitGuess), max(CoefLagFull)]+2)]);
            xlabel('X'); ylabel('\epsilon(x)'); grid on;   
            
            print('-depsc2',[folder,'/Run',num2str(RunIdx),'_Iter',num2str(n),'.eps']);

            
        end
    end
    Coef_local_1 = CoefLocal;   
    
end


% % propagate the data and then run the algorithm again: 
% if step <= 3
%     [DirDataNew, NeumanDataNew] = data_propagation_1d(CoefValues,SubGrid(1:NrRefine:end),t,DirichletData,NeumanData);
%     
%     
%     
% end



% run the local method alone: 
if step <= 4
    
    NrIter = 3;
    SubGrid = X_LagOld(2:end-1);
    InitGuess = (CoefLag1 - 1)*0; 

    ParBasis = FEM_basis(SubGrid,X_FDM); % the basis for parametrizing the coefficient
    CoefInitGuess = FEM_basis_expansion(ParBasis,InitGuess);


    for n = 1:NrIter
        fprintf('%s%d','Nr of refinements: ', n);
        lb = 0*InitGuess - 1; ub = lb + 30;

        [CoefLocal,CoefValues,fval,gradf] = waveeq1d_inverse_ls(Ui_tt,InitGuess,SubGrid,DirichletData,IdxMea,X_FDM,dt,lb,ub,options);

        % refine the mesh: uniform refinement
        SubGrid = linearinter(SubGrid,NrRefine);
        InitGuess = linearinter(CoefValues,NrRefine);

    %             % Test: refine the mesh at the points of large gradient only:
    %             SubGridOld = SubGrid; 
    %             SubGrid = mesh_refinement_1d(SubGridOld,gradf,0.85); 
    %             InitGuess = linearinterpolation(CoefValues,SubGridOld,SubGrid);            


        % --- display the coefficient: 
        figure; set(gca,'fontsize',16);
        hold off; plot(X_FDM,CoefLocal(idx),'--b','linewidth',2); 
        hold on; plot(X_FDM,CoefInitGuess(idx),'-k','linewidth',2); hold off; 
        hold on; plot(X_FDM,CoefExact(idx),'-.r','linewidth',2); hold off; 

        legend('Reconstruction using the local method','Initial guess','Exact coefficient');
        axis([X_FDM(1) X_FDM(end) min([min(CoefLocal), min(CoefExact), min(CoefInitGuess)])-0.5 max([max(CoefLocal), max(CoefExact), max(CoefInitGuess)]+2)]);
        xlabel('X'); ylabel('\epsilon(x)'); grid on;   

        print('-depsc2',[folder,'/local_method_only_iter',num2str(n),'.eps']);


    end
end


% if step == 5
%     disp('Running the local method: with Dirichlet and Neumann conditions');
%     % 2. Estimate the coefficient using the local method with adaptive parametrization
% 
%     % recursively adaptive parametrization:
%     NrIter = 1; % number of total refinements
%     NrRefine = 2; % refinement factor
%     NrUnknowns = length(CoefLag);
% 
%     % the first initial guess:
%     SubGrid = X_Lag(2:end-1);
%     InitGuess = CoefLag - 1; % start from the homogeneous medium
% 
% % %     ------------- for checking the accuracy of the gradient using the adjoint method:
% %     NrUnknowns = 2; 
% %     SubGrid = linspace(0.6, 0.7,NrUnknowns);
% %     InitGuess = 2*ones(NrUnknowns,1);
% % %     -------------
% 
%     ParBasis = FEM_basis(SubGrid,X_FDM); % the basis for parametrizing the coefficient
%     CoefInitGuess = FEM_basis_expansion(ParBasis,InitGuess);
% 
%     options = optimset('MaxIter',80,'TolFun',1e-6,'Algorithm','sqp'); % options for the fmin
% 
%     for n = 1:NrIter
%         fprintf('%s%d','Nr of refinements: ', n);
%         lb = 0*InitGuess; ub = lb + 7;
% 
%         % using both the Dirichlet and Neumann data for inversion: 
%         [Coef,CoefValues,fval] = waveeq1d_inverse_ls_DN_right(Ui_tt(:,1:IdxMea),InitGuess,SubGrid,DirichletData,NeumannData/2,X_FDM(1:IdxMea),dt,lb,ub,options);
%         Coef = [Coef; CoefExact(IdxMea+1:end)];
% 
%         
%         % refine the mesh:
%         SubGrid = linearinter(SubGrid,NrRefine);
%         InitGuess = linearinter(CoefValues,NrRefine);
% 
%         figure(n); set(gca,'fontsize',16);
%         hold off; plot(X_FDM,Coef,'--b','linewidth',2); 
%     %     hold on; plot(X_FDM,CoefInitGuess,'--r','linewidth',2); hold off; 
%         hold on; plot(X_FDM,CoefLagFull,'-k','linewidth',2); hold off; 
%         hold on; plot(X_FDM,CoefExact,'-.r','linewidth',2); hold off; 
% 
%         legend('Reconstructed using the local method','Reconstructed using Laguerre method','Exact coefficient');
%         axis([X_FDM(1) X_FDM(end) 0.5 max([max(Coef), max(CoefExact), max(CoefInitGuess), max(CoefLagFull)]+0.5)]);
%         xlabel('X'); ylabel('\epsilon(x)'); grid on;
%         title(n);
% 
%     end
% end




