function nzLocation = inferenceBlockHard(optIn)
    %parameters
    L = optIn.L;
    B = optIn.B;
    N = optIn.N;
    logRatio = optIn.logRatio;
    
    logRatio = max(-optIn.maxlogRatio, logRatio);
    logRatio = min(optIn.maxlogRatio, logRatio);
    
    %initialization
    nzLocation = zeros(N,L);
    pVector = [ones(1,B),zeros(1,L+B-2)];
    pMatrix = zeros(B + L - 1, L + 2 * B - 2);
    for i = 1:(L + B - 1)
        pMatrix(i,:) = circshift(pVector,i - 1);
    end
    
    for j = 1:N
        c = [zeros(B-1,1);logRatio(j,:)';zeros(B-1,1)];

       
        %normalization
        lc = pMatrix * c;
        lc0 = -max(lc);
        lc = lc - max(lc);
        plc = exp(lc);
        plc0 = exp(lc0);
        sumplc = sum(plc) + plc0;
        lc = log(plc ./ (sumplc - plc0));
        lc = max(-optIn.maxlogRatio, lc);
        lc = min(optIn.maxlogRatio, lc);
        
        loc = find(lc == max(lc));
        locEnd = min(loc(1), optIn.L);
        locStart = max(loc(1) - optIn.B + 1,1);
        nzLocation(j, locStart:locEnd) = 1; 
    end
end