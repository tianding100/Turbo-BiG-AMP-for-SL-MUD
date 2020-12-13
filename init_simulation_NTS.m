% Init for simulation
% Packets are QPSK modulated
% Generated by Poisson point process

clear;
clc;

setupPath;

setting.SNRList = -10:3:10;
numSNR = length(setting.SNRList);
setting.dataPath = 'data1.mat';
setting.randSeed = rng();
progress.idx = 1;   
progress.jdx = 1;

%% generate simulation parameters for nts MaDMA
%input parameters
optIn.M = 40;
optIn.L = 256;
optIn.winInt = 64;
optIn.B = 64;
optIn.U = 200;
optIn.Plambda = 2000;
optIn.Pnum = 100000;
optIn.alphabet = [1,-1,1i,-1i];
optIn.guardInterval = 64;
optIn.eps = 0.5;
optIn.spar = 0.2;
optIn.verbose = false;
optIn.maxTrials = 1;
optIn.p = 1;
optIn.blockInference = 1;
optIn.constellation = 1;
optIn.pathLossMax = 5; %the max difference of pathloss component in dB
setting.optIn = optIn;


% generate channel realizations and result file
genSimPar_NTS;
results = zeros(numPacket, numSNR);
progress.simPacket = zeros(caseNum,numSNR);
save(setting.dataPath,'setting','progress','results');