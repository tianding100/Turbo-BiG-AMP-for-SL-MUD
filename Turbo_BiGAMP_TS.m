function results = Turbo_BiGAMP_TS(optIn)
    if nargin == 0
%         load randseed_test.mat;
%         rng(s);
        
        optIn.maxTrials = 5;

      
        %Problem dimensions
        optIn.M = 40; %size of signal
        optIn.N = 40; %size of dictionary
        optIn.L = 256; %From Spielman et al.
        optIn.spar = 0.25;  % non-zero rate of each packet
        optIn.SNR = 15;
        optIn.alphabet = [1,-1,1i,-1i];
        optIn.verbose = true;
        
        optIn.rankContract = false;
        optIn.blockInference = 0;
              
        % Constellation type : 
        
        % 0 - Gaussian, 1 - QPSK
        optIn.constellation = 1;
        
        optIn.X = generateXiid(optIn);
        optIn.Avar = 1;
        optIn.A = (randn(optIn.M,optIn.N) + ...
            randn(optIn.M,optIn.N)*1i)/sqrt(2) * sqrt(optIn.Avar);
        optIn.W = (randn(optIn.M,optIn.L) + ...
            randn(optIn.M,optIn.L)*1i)/sqrt(2) * sqrt(10^(-optIn.SNR/10) * optIn.spar);
        optIn.Wvar = 10^(-optIn.SNR/10) * optIn.spar;
        
        
       %% Important:init
%         optIn.knowCSI = 1;
%         optIn.init = optIn.A;
%         
        %add large scale fading
%         delta = zeros(1,size(optIn.A,2));
%         %Generate Y and change N
%         for i = 1:size(optIn.A,2)
%             delta(i) = rand()*5 - 2.5;
%             delta(i) = 10^(delta(i)/10);
%             optIn.A(:,i) = optIn.A(:,i) * sqrt(delta(i));
%         end
%         
%         optIn.Avar = optIn.Avar * mean(delta);
        
    end

    %% Problem Setup
        %% Defualt Parameters
    optInDefault.maxTrials = 5;
    optInDefault.inIt = 200;
    optInDefault.outIt = 20;
    optInDefault.verbose = 0;
    optInDefault.constellation = 1;
    optInDefault.blockInference = 0;
    optInDefault.alphabet = [1,-1,1i,-1i];
    optInDefault.eps = 0.5;
    optInDefault.M = size(optIn.A,1);
    optInDefault.N = size(optIn.A,2);
    optInDefault.Nbar = [];
    optInDefault.knowCSI = 0;
    optInDefault.L = size(optIn.X,2);
    optInDefault.Y = optIn.A * optIn.X + optIn.W; 
    optInDefault.init = 'random';
    if ~isfield(optIn,'Wvar')
        optInDefault.Wvar = 10^(-optIn.SNR/10) * optIn.spar;
    end
    optInDefault.Avar = 1;
    optInDefault.rankContract = false;
    optIn = checkOptions(optInDefault,optIn);
    
    %Problem Setup
    maxTrials = optIn.maxTrials;
    M = optIn.M;
    L = optIn.L;
    N = optIn.N;
    spar = optIn.spar;
    X = optIn.X;
    A = optIn.A;
    W = optIn.W;
    Y = optIn.Y;
    
    Nbar = optIn.Nbar;
    if optIn.rankContract && isempty(Nbar)
        Nbar = ceil(mean(abs(Y(:)).^2)/optIn.Avar/spar * 1.2);
    else
        Nbar = N;
    end
    blockInference = optIn.blockInference;
    
    
    Y = optIn.Y;
    
    %set BiGAMP option
    opt = BiGAMPOpt();  %initialize the options object
    opt.nit = optIn.inIt;  %limit iterations
    opt.verbose = false;
    opt.blockInference = blockInference;
    
    
    problem = BiGAMPProblem();
    problem.M = M;
    problem.N = N;
    problem.L = L;
    problem.spar = spar;

    %set EMopt
    EMopt = [];
    EMopt.verbose = optIn.verbose;
    EMopt.noise_var = optIn.Wvar;
    EMopt.active_var = 1;
    EMopt.nuA = optIn.Avar;
    EMopt.lambda = spar;
    EMopt.maxEMiter = optIn.outIt;
    EMopt.init = optIn.init; %'data';
    EMopt.knowCSI = optIn.knowCSI;
    
    %%%note
    EMopt.constellation = optIn.constellation;
    EMopt.alphabet = optIn.alphabet;
    EMopt.rankContract = optIn.rankContract;
    

    %% BiGAMP
    
    results.PER = 1;
    SER = [];
    results.Xhat = zeros(size(X));
    AMSE = inf;
    results.recIdx = [];

    for trial = 1:maxTrials
        %Run BiGAMP  
        [estFinTemp,~,~] = ...
            EMBiGAMP(Y,problem,opt,EMopt);
        
        errRes = checkErrorBiGAMP_TS(optIn,estFinTemp);
        recIdx = errRes.recIdx;
        results.Xhat(errRes.recIdx,:) = errRes.Xhat(errRes.recIdx,:);
        AMSE = min(AMSE, errRes.AMSE);
        results.recIdx = unique([results.recIdx,recIdx]);
        SER = [SER, errRes.SER];
        results.PER = 1 - length(results.recIdx) / N;
        if results.PER == 0
            break;
        end
    end
    
   
    results.SER = mean(SER,2);
    results.SER(results.recIdx) = 0;
    results.AMSE = AMSE;

    %% Show Results

    if (~isfield(optIn,'verbose')) || (optIn.verbose)
        disp(results);
        subplot(2,1,1);
        drawMatrix(optIn.X);
        subplot(2,1,2);
        drawMatrix(results.Xhat);
    end
end