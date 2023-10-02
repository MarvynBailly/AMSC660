function result = cholesky(A)
    n = length(A);
    L = zeros(n,n);

    % Check if A is SPD
    for j = 1 : n
        L(j,j) = (A(j,j) - sum(L(j,1:j-1).^2))^(1/2); 

        if(L(j,j) == 0 || ~isreal(L(j,j)))
            fprintf("The matrix is not positive definite")
            result = 0;
            return 
        end
       
        for i = j + 1 : n
            L(i,j) = (A(i,j) - sum(L(i,1:j-1).*L(j,1:j-1)))/L(j,j); 
        end
    end
    result = L;
end
