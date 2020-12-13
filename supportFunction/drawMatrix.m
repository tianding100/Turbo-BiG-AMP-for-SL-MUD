function output = drawMatrix(X,XD)
    if (nargin == 1)
        A = 2 * abs(X);
        map = [1, 1, 1; 0, 0, 0];
        output = imagesc(A);
        colormap(map);
    else
        map = [1, 1, 1; 1, 0, 0 ; 0, 0, 0];
        %0--white, 1--red, 2--black
        errorFlag = (X~=XD) * 2;
        A = min(abs(X) + errorFlag, 2);
        A(A > 0) = 3 - A(A > 0);
        output = imagesc(A);
        colormap(map);
    end
    grid on;
end