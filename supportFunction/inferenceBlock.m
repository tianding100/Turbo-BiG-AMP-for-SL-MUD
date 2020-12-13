function newSparsity = inferenceBlock(optIn)
    
    %damping factor
    alpha = 1;

    %parameters
    L = optIn.L;
    B = optIn.B;
    N = optIn.N;
    logRatio = optIn.logRatio;
    
    logRatio = max(-optIn.maxlogRatio, logRatio);
    logRatio = min(optIn.maxlogRatio, logRatio);
    
    %initialization
    lx = zeros(N,L);
    pVector = [ones(1,B),zeros(1,L+B-2)];
    pMatrix = zeros(B + L - 1, L + 2 * B - 2);
    
    for i = 1:(L + B - 1)
        pMatrix(i + 1,:) = circshift(pVector,i - 1);
    end
    
    Ctrunc = 1 - eye(L);
    
    for j = 1:N   
      C = [zeros(B-1,L);repmat(logRatio(j,:)',1,L).*Ctrunc;zeros(B-1,L)];
      LC = pMatrix * C;
      LC = max(-optIn.maxlogRatio, LC);
      LC = min(optIn.maxlogRatio, LC);
      S1 = sum(exp(LC) .* pMatrix(:,B:(B+L-1)));
      S0 = sum(exp(LC) .* (1 - pMatrix(:,B:(B+L-1))));
      lx(j,:) = log(S1) - log(S0);
    end
        
    
    lx = alpha * (lx - logRatio) + logRatio;
    
    lx = max(-optIn.maxlogRatio, lx);
    lx = min(optIn.maxlogRatio, lx);
    
    newSparsity = 1 ./ (1 + exp(-lx));
end