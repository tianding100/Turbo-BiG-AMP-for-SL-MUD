function results = Turbo_BiGAMP_NTS(optIn)
    setupPath;
    if nargin == 0
        %example problem setup  
        randSeed = rng;
        save rs.mat randSeed
%         load rs.mat;
%         rng(randSeed);
        
        optIn.maxTrials = 1;
        
        optIn.M = 40; %  size of signal
        optIn.N = 40; %  size of dictionary
        optIn.Nbar = 45; % (estimated) number of active users
        optIn.B = 64;
        optIn.L = 256; 
        optIn.spar = 0.2;  % sparsity level of each packet
        optIn.SNR = 2;
        optIn.alphabet = [1,-1,1i,-1i];
        optIn.verbose = true;
        optIn.rankContract = true;
        
        % Constellation type : 
        % 0 - Gaussian, 1 - QPSK
        optIn.constellation = 1;

        % blockInference on/off
        optIn.blockInference = 1;
        
        optIn.verbose = true;
        
        %generate data sample
        optIn.X = generateX(optIn);
        optIn.Avar = 1;
        optIn.A = (randn(optIn.M,optIn.N) + ...
            randn(optIn.M,optIn.N)*1i)/sqrt(2) * sqrt(optIn.Avar);
        optIn.W = (randn(optIn.M,optIn.L) + ...
            randn(optIn.M,optIn.L)*1i)/sqrt(2) * sqrt(10^(-optIn.SNR/10));
        optIn.Wvar = 10^(-optIn.SNR/10);
        
        %add large scale fading
        delta = zeros(1,size(optIn.A,2));
        %Generate Y and change N
        for i = 1:size(optIn.A,2)
            delta(i) = rand()*10 - 5;
            delta(i) = 10^(delta(i)/10);
            optIn.A(:,i) = optIn.A(:,i) * sqrt(delta(i));
        end
        
        optIn.Avar = optIn.Avar * mean(delta);
        
        optIn.Y = optIn.A * optIn.X + optIn.W;
        
       %% Important:init
%         optIn.knowCSI = 1;
%         optIn.init = optIn.A;
    end
    
    %% Defualt Parameters
    optInDefault.maxTrials = 5;
    optInDefault.inIt = 200;
    optInDefault.outIt = 20;
    optInDefault.verbose = 0;
    optInDefault.constellation = 1;
    optInDefault.blockInference = 1;
    optInDefault.alphabet = [1,-1,1i,-1i]';
    optInDefault.eps = 0.5;
    optInDefault.M = size(optIn.A,1);
    optInDefault.N = size(optIn.A,2);
    optInDefault.Nbar = [];
    optInDefault.L = size(optIn.X,2);
    optInDefault.Y = optIn.A * optIn.X + optIn.W; 
    if ~isfield(optIn,'Wvar')
        optInDefault.Wvar = 10^(-optIn.SNR/10);
    end
    optInDefault.Avar = 1;
    optInDefault.rankContract = true;
    optInDefault.init = 'random';
    optInDefault.knowCSI = 0;
    optIn = checkOptions(optInDefault,optIn);
    
    %% Problem Setup
    maxTrials = optIn.maxTrials;
    Y = optIn.Y;
    
    %estimated no. users
    Nbar = optIn.Nbar;
    if optIn.rankContract && isempty(Nbar)
        Nbar = ceil(mean(abs(Y(:)).^2)/optIn.Avar/optIn.spar * 1.2);
    end
    
    
    %set BiGAMP option
    opt = BiGAMPOpt();  %initialize the options object
    opt.nit = optIn.inIt;  %limit iterations
    opt.verbose = false;
    opt.blockInference = optIn.blockInference;
    
    
    problem = BiGAMPProblem();
    problem.M = optIn.M;
    problem.N = Nbar;  %set initial rank
    problem.L = optIn.L;
    problem.B = optIn.B;
    problem.spar = optIn.spar;
    
    
    %set EMopt
    EMopt = [];
    EMopt.verbose = optIn.verbose;
    EMopt.noise_var = optIn.Wvar;
    EMopt.active_var = 1;
    EMopt.nuA = optIn.Avar;
    EMopt.lambda = optIn.spar;
    EMopt.maxEMiter = optIn.outIt;
    EMopt.constellation = optIn.constellation;
    EMopt.alphabet = optIn.alphabet;
    EMopt.rankContract = optIn.rankContract;
    EMopt.init = optIn.init; %'data';
    EMopt.knowCSI = optIn.knowCSI;
    
    %% BiGAMP
    results.PER = 1;
    results.MSE = [];
    results.recIdx = [];
    results.numPac = size(optIn.X,1);
    results.Xhat = zeros(size(optIn.X));

    for trial = 1:maxTrials
        %Run BiGAMP  
        [estFinTemp,~,~] = ...
            EMBiGAMP(Y,problem,opt,EMopt);
        
        errRes = checkErrorBiGAMP_NTS(optIn,estFinTemp);
        recIdx = errRes.recIdx;
        results.numPac = errRes.numPac;
        results.Xhat = errRes.Xhat;
        results.MSE = errRes.MSE;
        results.recIdx = unique([results.recIdx,recIdx]);
        results.PER = 1 - length(results.recIdx) / max(results.numPac,1);
        results.Nbar = size(estFinTemp.xhat,1);
        if results.PER == 0
            break;
        end
    end
    
    %% Show Results

    if (~isfield(optIn,'verbose')) || (optIn.verbose)
        disp(results);
        subplot(2,1,1);
        drawMatrix(optIn.X);
        subplot(2,1,2);
        drawMatrix(results.Xhat);
    end
end