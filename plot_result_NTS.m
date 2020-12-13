%plot simulation results
setupPath;
dataPath = 'data1.mat';

%calculate PER
load(dataPath);
rng(setting.randSeed);
optIn = setting.optIn;
genSimPar_NTS;
pNum = Signal.winTable(3,progress.idx - 1) - 1;
y = 1 - (sum(results(1:pNum,:))/pNum);
x = setting.SNRList;

%plot
h1 = semilogy(x,y,'--','Marker','o','Color','blue','MarkerSize',12,'linewidth',1);
grid on;

ylabel('Packet error rates');
xlabel('P (dBm)');
set(gcf,'color','white');
set(gca,'FontSize',14);
grid on;