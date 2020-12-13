%generate signals, channels and noise
Signal = generateSignal_NTS(optIn);
symTime = Signal.timeTable(end);    %timetable of user packets
numPacket = length(Signal.timeTable);   %total number of user packets
caseNum = size(Signal.winTable,2);  
optIn.Avar = Signal.Avar;

Channel = (randn(optIn.M,numPacket)+randn(optIn.M,numPacket)*1i) / sqrt(2);
Channel = sqrt(Signal.betaTable).*Channel;
Noise = (randn(optIn.M,symTime + optIn.L)+randn(optIn.M,symTime + optIn.L)*1i) / sqrt(2) ;
simOptIn = optIn;
numPacket = length(Signal.timeTable);   %total number of user packets