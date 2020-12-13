function errRes = checkErrorBiGAMP_TS(optIn,est)
    %recover signals from optIn
    A = optIn.A;
    X = optIn.X;
    N = optIn.N;
    L = optIn.L;
    eps = optIn.eps;
    alphabet = optIn.alphabet;
    
    %1 - Determine zero entries
    Xhat0 = est.xhat;
    Ahat0 = est.Ahat;
    Nhat = size(Xhat0,1);
    Xhat1 = Xhat0 .* (abs(Xhat0) > eps);
    
    %2 - Eliminate phase ambiguousness
    for i = 1:Nhat
        temp = find(abs(Xhat1(i,:)) > 0);
        if ~isempty(temp)
            %First symbol is 1
            scale = Xhat1(i,temp(1)) / abs(Xhat1(i,temp(1)));
            Xhat1(i,:) = Xhat1(i,:) ./ scale;
            Ahat0(:,i) = Ahat0(:,i) .* scale; 
            Xhat0(i,:) = Xhat0(i,:) ./ scale;
        end
    end
    
    %3 - Hard decision
    for i = 1:Nhat
        for j = 1:L
            if abs(Xhat0(i,j)) > 0
                [~,loc] = min(abs(alphabet - Xhat0(i,j)));
                Xhat0(i,j) = alphabet(loc);
            end
        end
    end
        
    %Record the recovered packets
    Xhat = zeros(size(X));
    Ahat = zeros(size(A));
    recFlag = zeros(1,N);
    idx = 1:1:Nhat;
    for i = 1:N
        %find the row in Xhat0 closest to X(i,:)
        sumErr = L;
        tempid = 1;
        for j = 1:length(idx)
            newErr = sum(Xhat0(idx(j),:) ~= X(i,:));
            if newErr < sumErr
                tempid = j;
                sumErr = newErr;
                if (sumErr == 0)
                    recFlag(i) = idx(j);
                    break;
                end
            end
        end
        if ~isempty(idx)
            Xhat(i,:) = Xhat0(idx(tempid),:);
            Ahat(:,i) = Ahat0(:,idx(tempid));
            idx(tempid) = [];
        end
    end    
    
    
    errRes.recIdx = find(recFlag > 0);
%    Q = find_permutation(A,Ahat);
    %||A - BQ||
%    Ahat = Ahat * Q;
    errRes.AMSE = norm(Ahat - A,'fro')^2/size(A,1)/size(A,2);
    errRes.SER = sum(Xhat~=X,2)/ L;
    errRes.Xhat = Xhat;
    errRes.Ahat = Ahat;
    errRes.PER = 1 - length(errRes.recIdx)/N;
end