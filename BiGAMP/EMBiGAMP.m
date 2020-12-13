function [estFin,optFin,EMoptFin,estHist]...
    = EMBiGAMP(Y,problem,BiGAMPopt, EMopt)

% INPUTS:
% -------
% Y: the noisy data matrix
% problem: An objet of the class BiGAMPProblem specifying the problem
%   setup, including the matrix dimensions and observation locations
% opt (optional):  A set of options of the class BiGAMPOpt.
% EMopt (optional): A structure containing several fields. If the input is
%   omitted or empty, the values in defaultOptions defined below will be
%   used. The user may pass a structure containing any or all of these fields
%   to replace the defaults. Fields not included in the user input will be
%   replaced with the defaults
%
% OUTPUTS:
% --------
% estFin: Structure containing final BiG-AMP outputs
% optFin: The BiGAMPOpt object used
% EMoptFin: The EM options used
% estHist: Structure containing per iteration metrics about the run

%% Handle Options

%Create options if not provided
if (nargin < 3)
    BiGAMPopt = [];
end

%Handle empty options object
if (isempty(BiGAMPopt))
    BiGAMPopt = BiGAMPOpt();
end

%If no inputs provided, set user options to empty
if nargin < 4
    EMopt = [];
end


%Get problem dimensions
M = problem.M;
L = problem.L;
N = problem.N;

if isfield(EMopt,'nuA') && ~isempty(EMopt.nuA)
    nuA = EMopt.nuA;
end


%Initial values (note that the mean of the X entries is initialized
%with zero)
defaultOptions.noise_var = []; %initial noise variance estimate
defaultOptions.active_var = []; %initial variance estimate of active X elements

%Options to control learning
defaultOptions.learn_noisevar= false; %learn the variance of the AWGN
defaultOptions.learn_var = false; %learn the variance of the active X entries
defaultOptions.learn_mean = false; %learn the mean of the active X entries
defaultOptions.learn_lambda = false; %learn the sparsity of X entries
defaultOptions.sig_dim = 'joint'; %learn a single variances for X entries (joint)
%or a different variance per (row) or (column)

%Iteration control
defaultOptions.maxEMiter = 20; %maximum number of EM cycles
defaultOptions.tmax = 20; %first EM iteration to use full expression for noise variance update

%Tolerances
defaultOptions.EMtol = BiGAMPopt.tol; %convergence tolerance
defaultOptions.maxTol = 1e-4;%largest allowed tolerance for a single EM iteration

%Init method for dictionary (choices are 'data', 'random', and 'SVD')
defaultOptions.init = 'random';

%Warm starting, set to true to enable. Warm starting may perform
%faster, but somewhat increases the likelihood of getting stuck
defaultOptions.warm_start = false;

%rank selection parameters
defaultOptions.tau_mos = 5;
rankSelFlag = 0;

%Combine user options with defaults. Fields not provided will be given
%default values
EMopt = checkOptions(defaultOptions,EMopt);


%% Initial Setup

%Indices of observed entries- subsampling has not been tested with DL
if ~isempty(problem.rowLocations)
    omega = sub2ind([M L],problem.rowLocations,problem.columnLocations);
else
    omega = true(M,L);
end

%Sparse mode not yet allowed
if BiGAMPopt.sparseMode
    error('Sparse mode not yet supported');
end


%Set noise variance if requested
if ~isfield(EMopt,'noise_var') || isempty(EMopt.noise_var)
    if ~BiGAMPopt.sparseMode
        EMopt.noise_var = sum(sum(abs(Y(omega)).^2))/(sum(sum(omega))*101);
    else
        EMopt.noise_var = sum(sum(abs(Y).^2))/...
            (M*L*101);
    end
end


%History
histFlag = false;
if nargout >=4
    histFlag = true;
    estHist.errZ = [];
    estHist.errX = [];
    estHist.errA = [];
    estHist.val = [];
    estHist.step = [];
    estHist.pass = [];
    estHist.timing = [];
end

%Set initial noise variance
nuw = EMopt.noise_var;

%Init X params
meanX = zeros(N,L);

if ~isfield(EMopt,'lambda') || isempty(EMopt.lambda)
    lambda = 0.1*ones(N,L);
else
    lambda = EMopt.lambda;
end


%Compute X variance, setting nuA = 1
if ~isfield(EMopt,'active_var') || isempty(EMopt.active_var)
    nuX = (norm(Y,'fro')^2/numel(omega) - nuw)/N/mean(lambda(:));
else
    nuX = EMopt.active_var;
end


%Specify tmax
if ~isfield(EMopt,'tmax') || isempty(EMopt.tmax)
    tmax = EMopt.maxEMiter;
else
    tmax = EMopt.tmax;
end


%Initialize loop
t = 0;
stop = 0;
state = [];

%Initialize xhat with zeros
xhat = zeros(N,L); %seems to perform better


%Initialize the dictionary
if isnumeric(EMopt.init)
    Ahat = EMopt.init;
    W_LMMSE = Ahat' * pinv(Ahat * Ahat' + EMopt.noise_var * size(length(Y),1));
    xhat = W_LMMSE * Y;
else
    if strcmp(EMopt.init,'data') %init from the data
        whichCol = 0;
        drawAttempts = 0;
        Ycond = cond(Y);
        while(whichCol <= N) && (drawAttempts < 20)
            
            %Init
            atomOrder = randperm(L); %different random order each time
            Ahat = zeros(M,N);
            Ahat(:,1) = Y(:,atomOrder(1))/norm(Y(:,atomOrder(1)));
            whichCol = 2;
            counter = 2;
            
            %Try to draw the initial dictionary
            while whichCol <= N && counter <= L
                
                %Assign the new column
                Ahat(:,whichCol) = Y(:,atomOrder(counter)) /...
                    norm(Y(:,atomOrder(counter)));
                
                %Check inner product and increment if different
                %also check cond to watch for rank deficient cases
                if (max(abs(Ahat(:,1:(whichCol-1))'*Ahat(:,whichCol))) < 0.9) && ...
                        (cond(Ahat(:,1:whichCol)) < 10*Ycond)
                    whichCol = whichCol + 1;
                end
                
                %Increment counter
                counter = counter + 1;
            end
            
            drawAttempts = drawAttempts + 1;
            
        end
        if (EMopt.verbose)
         disp(['Dictionary draw finished after attempts: '...
             num2str(drawAttempts)]);
        end
%         
        %Fill in with random if needed
        if whichCol <= N
            
            Ahat(:,whichCol:end) = randn(size(Ahat(:,whichCol:end))).*sqrt(nuA);
        end
        
    elseif strcmp(EMopt.init,'random') %or random
        
        %Draw random
        Ahat = (randn(M,N) + randn(M,N)*1i) / sqrt(2) .* sqrt(nuA);
        
    elseif strcmp(EMopt.init,'SVD') %or SVD
        
        %Use the econ SVD
        [U_1,S_1,V_1] = svd(Y,'econ');
        Ahat = zeros(M,N);
        Ahat(:,1:M) = U_1*S_1;
        
        %Scale Ahat so that it has unit variance entries
        vval = var(reshape(Ahat(:,1:M),1,[]));
        Ahat = Ahat / sqrt(vval);
        
        %Fill the rest of A with random entries
        Ahat(:,(M+1):end) = randn(M,(N - M));
        
        %Assign xhat to make the SVD exact
        xhat = zeros(N,L);
        xhat(1:M,:) = V_1'*sqrt(vval); 
    end
end

%Set init
BiGAMPopt.xhat0 = xhat;
BiGAMPopt.Ahat0 = Ahat;
BiGAMPopt.Avar0 = ones(M,N) .* nuA;
BiGAMPopt.xvar0 = ones(N,L) .* nuX;

%Ensure that diagnostics are off to save run time
BiGAMPopt.diagnostics = 0;

%Ensure that EM outputs are calculated
BiGAMPopt.saveEM = 1;

%Original tolerance
tol0 = BiGAMPopt.tol;


%Outer loop for rank learning
SNR = 100;
zhatOld = 0;


%% Main Loop

%EM iterations
%The < 2 condition is so that we can do a final iteration using the full
%noise varaince update without repeating code
while stop < 2
    
    %Start timing
    tstart = tic;
    
    %Estimate SNR
    SNRold = SNR;
    if t > 0
        SNR = norm(zhat,'fro')^2 / norm(Y - zhat,'fro')^2;
    end
    %Set tolerance
    %tolNew = 1 / SNR;
    tolNew = min(max(tol0,1/SNR),EMopt.maxTol);
    BiGAMPopt.tol = tolNew;
    
    %Increment time exit loop if exceeds maximum time
    t = t + 1;
    if t >= EMopt.maxEMiter || stop > 0
        stop = stop + 1;
    end
    
    %Prior on A
    if (EMopt.knowCSI)
        gA = CFixEstimIn(Ahat);
    else
        gA = CAwgnEstimIn(0, EMopt.nuA);
    end
    
    %Prior on X
    if (EMopt.constellation == 1) && (BiGAMPopt.blockInference)
        gXbase = QPSKEstimIn_NTS();
    elseif (EMopt.constellation == 1)
        gXbase = QPSKEstimIn_TS();
    else
        gXbase = CAwgnEstimIn(meanX, nuX);
    end
    gX = SparseEstim(gXbase,lambda);

    %Output log likelihood
    gOut = CAwgnEstimOut(Y, nuw);
    
    %Stop timing
    t1 = toc(tstart);
    
    %Clear the state if not warm starting
    if ~EMopt.warm_start
        state = [];
    end
    
    %Run BiG-AMP
    [estFin2,~,estHist2,state] = ...
        BiGAMP(gX, gA, gOut, problem, BiGAMPopt,state);
    
    %Start timing
    tstart = tic;
    
    %Correct cost function
    estHist2.val = estHist2.val -0.5*numel(omega)*log(2*pi*nuw);
    
    
    %Determine sparsity pattern
    [~,~,~,p1] = gX.estim(estFin2.rhat,estFin2.rvar);
    p1 = abs(p1);
    
    %Report progress
    if histFlag
        error_value = estHist2.errZ(end);
    else
        error_value = nan;
    end
    if EMopt.verbose
         disp(['It ' num2str(t,'%04d')...
             ' nuX = ' num2str(mean(nuX(:)),'%5.3e')...
            ' Lam = ' num2str(mean(lambda(:)),'%0.2f')...
            ' tol = ' num2str(BiGAMPopt.tol,'%5.3e')...
            ' SNR = ' num2str(10*log10(SNR),'%03.2f')...
            ' Z_e = ' num2str(error_value,'%05.4f')...
            ' nuw = ' num2str(nuw,'%5.3e')...
            ' Avg Spar = ' num2str(sum(sum(p1))/L,'%3.1f')...
            ' numIt = ' num2str(length(estHist2.errZ),'%04d')])
    end    
    %Compute zhat
    zhat = estFin2.Ahat*estFin2.xhat;
    
    %Calculate the change in signal estimates
    norm_change = norm(zhat-zhatOld,'fro')^2/norm(zhat,'fro')^2;
    
    %Check for estimate tolerance threshold
    if (norm_change < max(tolNew/10,EMopt.EMtol)) &&...
            ( (norm_change < EMopt.EMtol) ||...
            (abs(10*log10(SNRold) - 10*log10(SNR)) < 1))
        stop = stop + 1;
    end
    
    %Update noise variance. Include only a portion of the Zvar
    %in beginning stages of EMGMAMP because true update may make it
    %unstable.
    %Learn noise variance
%     if EMopt.learn_noisevar
%         
%         %First, just the term based on the residual
%         if BiGAMPopt.sparseMode
%             nuw = norm(Y - zhat,'fro')^2 ...
%                 /numel(omega);
%         else
%             
%             nuw = norm(Y(omega) - zhat(omega),'fro')^2 ...
%                 /numel(omega);
%         end
%         
%         %Then the component based on zvar
%         if t >= tmax || stop > 0
%             if isscalar(estFin2.zvar)
%                 nuw = nuw + estFin2.zvar;
%             else
%                 if BiGAMPopt.sparseMode
%                     nuw = nuw + sum(estFin2.zvar)/numel(omega);
%                 else
%                     nuw = nuw + sum(estFin2.zvar(omega))/numel(omega);
%                 end
%             end
%         end
%     end
    
    %Estimate new X parameters by EM
%     [lambda, meanX, nuX] =...
%         BG_update(estFin2.rhat, estFin2.rvar,...
%         lambda, meanX, nuX, EMopt);
    
    %Reinitialize GAMP estimates
    zhatOld = zhat;
    BiGAMPopt.Ahat0 = estFin2.Ahat;
    BiGAMPopt.Avar0 = estFin2.Avar;
    BiGAMPopt.step = BiGAMPopt.stepMin; %estHist2.step(end);
    
    %rank contraction
    if EMopt.rankContract
        Flagmax = 3;
        if (rankSelFlag < Flagmax)
            singval = svd(estFin2.xhat);
            if singval(1) > 0.001
                N = length(singval);
                Rval = singval(1:(N-1)) ./ singval(2:N);
                Nnew = find(Rval == max(Rval));
                if (~isempty(Nnew))
                    Nnew = Nnew(1);
                    if Rval(Nnew) > (sum(Rval)-Rval(Nnew)) * EMopt.tau_mos/(N-2)
                        %select the new rank
                        rankSelFlag = rankSelFlag + 1;
                        %prune matrices
                        problem.N = Nnew;
                        BiGAMPopt.xhat0 = BiGAMPopt.xhat0(1:Nnew,:);
                        BiGAMPopt.xvar0 = BiGAMPopt.xvar0(1:Nnew,:);
                        BiGAMPopt.Ahat0 = BiGAMPopt.Ahat0(:,1:Nnew);
                        BiGAMPopt.Avar0 = estFin2.Avar(:,1:Nnew);
                        BiGAMPopt.step = BiGAMPopt.stepMin;
                    end
                end
            end
        end
    end    
    %Stop timing
    t2 = toc(tstart);
    
    %Output Histories if necessary
    if histFlag
        estHist.errZ = [estHist.errZ; estHist2.errZ];
        estHist.errX = [estHist.errX; estHist2.errX];
        estHist.errA = [estHist.errA; estHist2.errA];
        estHist.val = [estHist.val; estHist2.val];
        estHist.step = [estHist.step; estHist2.step];
        estHist.pass = [estHist.pass; estHist2.pass];
        if t == 1
            estHist.timing = [estHist.timing; t1 + t2 + estHist2.timing];
        else
            estHist.timing = [estHist.timing; t1 + t2 + estHist.timing(end) ...
                + estHist2.timing];
        end
    end
end

%% Cleanup

%Update finals
estFin = estFin2;
optFin = BiGAMPopt;
EMoptFin = EMopt;

if BiGAMPopt.blockInference
    estFin.xhat = estFin.xhat .* estFin.nzLocation;
end



%Include learned parameters
if histFlag
    estHist.meanX = meanX;
    estHist.nuX = nuX;
    estHist.lambda = lambda;
    estHist.nuw = nuw;
    estHist.p1 = p1;
end




