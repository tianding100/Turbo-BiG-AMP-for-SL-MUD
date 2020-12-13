clear;
clc;
setupPath;
dataPath = 'data1.mat';

%%
%load settings
load(dataPath);
SNRList = setting.SNRList;
numSNR = length(SNRList);
rng(setting.randSeed);
optIn = setting.optIn;
genSimPar_NTS;
optIn = simOptIn;
warning off;

%%
while progress.idx <= caseNum
    idx = progress.idx;
    X = strcat('case',num2str(idx));
    disp(X);
    
    t0 = Signal.winTable(1,idx);
    t1 = Signal.winTable(2,idx);
    flag0 = Signal.winTable(3,idx);
    flag1 = Signal.winTable(4,idx);
    
    %construct the problem model
    optIn.X = zeros(flag1 - flag0 + 1, optIn.L);
    optIn.A = zeros(optIn.M,flag1 - flag0 + 1);
    for ii = flag0:flag1
        xStart = max(1,Signal.timeTable(ii) - t0 + 1);
        xEnd = min(optIn.L,Signal.timeTable(ii) - t0 + optIn.B);
        pStart = max(1,t0 - Signal.timeTable(ii) + 1);
        pEnd = pStart + xEnd - xStart;
        optIn.X(ii - flag0 + 1,xStart:xEnd) = Signal.packets(ii,pStart:pEnd);
        optIn.A(1:optIn.M,ii - flag0 + 1) = Channel(1:optIn.M,ii) * sqrt(Signal.betaTable(ii));
    end
    optIn.N = flag1 - flag0 + 1;
    PER = zeros(numSNR,1);
    randSeed = rng();
    while (progress.jdx <= numSNR)
        rng(randSeed);
        jdx = progress.jdx;

        optIn.SNR = SNRList(jdx);
        optIn.W = Noise(1:optIn.M,t0:t1) * sqrt(10^(-SNRList(jdx)/10));
        optIn.Wvar = 10^(-SNRList(jdx)/10);

        progress.simPacket(flag0:flag1,jdx) = 1;
        tempRes = Turbo_BiGAMP_NTS(optIn);
        windowPER = tempRes.PER;
        recIndex = tempRes.recIdx + flag0 - 1;
        results(recIndex,jdx) = 1;

        % display and save results      
        X = ['SNR:',num2str(SNRList(jdx)),'dB',' --- numPac:',num2str(flag0)];
        PER(jdx) = 1 - sum(results(1:flag0,jdx))/sum(progress.simPacket(1:flag0,jdx));
        X = [X,'--- PER:',num2str(PER(jdx))];
        disp(X);        
        progress.jdx = progress.jdx + 1;
        %save(dataPath,'setting','progress','results');
    end
    progress.jdx = 1;
    progress.idx = progress.idx +  1;
    save(dataPath,'setting','progress','results');
end



