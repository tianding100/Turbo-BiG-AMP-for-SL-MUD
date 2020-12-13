function errRes = checkErrorBiGAMP_NTS(optIn,est)
    %validate the recovered type-I packet and compute PER
    X = optIn.X;
    A = optIn.A;
    N = optIn.N;
    L = optIn.L;
    B = optIn.B;
    eps = optIn.eps;
    alphabet = optIn.alphabet;
    
    %1 - Determine zero entries
    Xhat0 = est.xhat;
    Ahat = est.Ahat;
    Nhat = size(Xhat0,1);
    Xhat1 = Xhat0 .* (abs(Xhat0) > eps);
    
    %2 - Eliminate phase ambiguousness
    for i = 1:Nhat
        temp = find(abs(Xhat1(i,:)) > 0);
        if ~isempty(temp)
            %First symbol is 1
            scale = Xhat1(i,temp(1)) / abs(Xhat1(i,temp(1)));
            Xhat1(i,:) = Xhat1(i,:) ./ scale;
            Ahat(:,i) = Ahat(:,i) .* scale; 
        end
    end

    %3 - Hard decision
    for i = 1:Nhat
        for j = 1:L
            if abs(Xhat1(i,j)) > 0
                [temp,loc] = min(abs(alphabet - Xhat1(i,j)));
                Xhat1(i,j) = alphabet(loc);
            end
        end
    end
    
    %4 - Compare the type-I packets
    %Extract the original type-I packets
    pacIdx = [];
    for i = 1:N
        if sum(X(i,:) ~= 0) == B
            pacIdx = [pacIdx;i];
        end
    end
    
    %Record the recovered type-I packets
    recFlag = zeros(1,N);
    for i = 1:Nhat
        %find the type-I packet in X closest to XhatHard(i,:)
        if (sum(Xhat1(i,:) ~= 0) == B)
            for j = 1:length(pacIdx)
                if sum(Xhat1(i,:) ~= X(pacIdx(j),:)) == 0
                   recFlag(pacIdx(j)) = i;
                   break;   
                end
            end
        end
    end
    
    Xhat = zeros(size(X));
    MSE = [];
    recIdx = find(recFlag > 0);
    for i = 1:length(recIdx)
        Xhat(recIdx(i),:) = Xhat1(recFlag(recIdx(i)),:);
        MSETemp = abs(Ahat(:,recFlag(recIdx(i))) - A(:,recIdx(i))).^2;
        MSE = [MSE;mean(MSETemp)];
    end
    
    errRes.numPac = length(pacIdx);
    errRes.Xhat = Xhat;
    errRes.MSE = MSE;
    errRes.recIdx = recIdx;
    errRes.PER = 1 - sum(recFlag ~= 0)/max(errRes.numPac,1);
end