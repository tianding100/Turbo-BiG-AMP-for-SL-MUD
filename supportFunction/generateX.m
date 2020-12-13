function X  = generateX(optIn)
    %generate X with consecutive locations of non-zero elements
    
    %parameters
    N = optIn.N;
    L = optIn.L;
    B = optIn.B;
    alphabet = optIn.alphabet;
    alLen = size(alphabet,2);
    loc = zeros(N,1);
    
    %determin locations of non-zero elements
    X = zeros(N,L);
    for i = 1:N
        loc(i) = ceil(rand() * (L + B)) - B;
        if (loc(i) > (1-B))
            X(i,max(loc(i),1):min(L,loc(i) + B - 1)) = 1;
        end
    end
    
    %generate modulated signals i.i.d.
    randMod = ceil(rand(N,L).*alLen);
    X = X.* alphabet(randMod);
    
    %The first bit of packet is 1
    for i = 1:N
        if (loc(i) >= 1)
            X(i,loc(i)) = 1;
        end
    end
end
