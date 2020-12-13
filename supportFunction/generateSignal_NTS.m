function Signal  = generateSignal_NTS(optIn)
    %generate X with consecutive locations of non-zero elements
    
    %parameters
    U = optIn.U;    %number of users
    L = optIn.L;        %window size
    B = optIn.B;        %packet size
    Plambda = optIn.Plambda;   %poisson distribution parameter
    Pnum = optIn.Pnum;      %total_number of package required   
    alphabet = optIn.alphabet;  %alphabet size
    guardInterval = optIn.guardInterval;    %guard interval between packets
    alLen = size(alphabet,2);      %constellation size
    winInt = optIn.winInt;
    
    uBeta = ones(U,1);
    for i = 1:U
        uBeta(i) = rand() * optIn.pathLossMax - optIn.pathLossMax/2;
        uBeta(i) = 10^(uBeta(i)/10);
    end
    Signal.Avar = mean(uBeta);    %no large fading
    
    %uTime first col is the time interval
    %uTime second col is the corresponding beta
    uTime = [zeros(U,1),uBeta];
    uTime(:,1) = ceil(exprnd(Plambda , U , 1));
    uTime(:,1) = uTime(:,1) - min(uTime(:,1)) + 1;   %next packet of each user
    uTime(:,1) = sort(uTime(:,1),1);
    
    PCount = 0;     %number of packets generated so far
    t0 = 1;     %start time of the considered window
    flag0 = 1;
    timeTable = [];
    betaTable = [];
    while PCount < Pnum
       t1 = uTime(1,1);
       beta = uTime(1,2);
       timeTable = [timeTable,t1];
       betaTable = [betaTable,beta];

       %update window time
       while (t1 - t0 >= L)
           flag0 = flag0 + 1;
           t0 = timeTable(flag0);
       end
       
       uTime(1,:) = [];
       uInt = ceil(exprnd(Plambda , 1 , 1));
       while (uInt < guardInterval)
           uInt = ceil(exprnd(Plambda , 1 , 1));
       end
       
       t1 = t1 + uInt;
       flag1 = find(uTime(:,1) > t1);
       if isempty(flag1)
           flag1 = size(uTime,1)+1;
       else
           flag1 = flag1(1);
       end
       newTimeRow = [t1,beta];
       uTime = [uTime(1:(flag1 - 1),:);newTimeRow;uTime(flag1:end,:)];
       PCount = length(timeTable);
    end
    
    
    %generate timetable for each window
    flag0 = 1;
    flag1 = 1;
    winTable = [];
    t0 = 1;
    t1 = L;
    symTime = timeTable(end);
    while t0 + L < symTime + B
        stuckFlag = 1;
        while stuckFlag
            flag0 = find(timeTable > t0 - B);
            flag0 = flag0(1);
            flag1 = find(timeTable <= t1);
            flag1 = flag1(end);
            if flag0 > flag1
                %no packets from t0 - t1
                t0 = timeTable(flag0);
                t1 = t0 + L - 1;
            else
                stuckFlag = 0;
            end
        end
        temp = [t0;t1;flag0;flag1];
        winTable = [winTable,temp];
        t0 = t0 + winInt;
        t1 = t1 + winInt;
    end
    
    Signal.timeTable = timeTable;
    Signal.packets = ceil(rand(length(timeTable),B) * alLen);
    Signal.packets = alphabet(Signal.packets);
    Signal.packets(:,1) = 1;
    Signal.winTable = winTable;
    Signal.betaTable = betaTable;
end
