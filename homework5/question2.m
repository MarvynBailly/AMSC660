function question2()
    At = rand(100);  %generate random guy

    A = At + At';    %compute symmetric matrix
    %A = At'*At;      %compute SPD matrix

    L = cholesky(A);
    if(L ~= 0)
        fprintf("The matrix is postive definte\n")
        mineval = min(eig(A));  %get smallest eigevalue
        if(mineval ~= 0 && isreal(mineval))
            fprintf("The minimal eigenvalue is %d which is positive and real\n",mineval)
        else
            fprintf("The algorithm failed")
        end 
        err = norm(L - chol(A));
        fprintf("The norm error of the algorithm and matlab cholesky is %d\n",err)
    end
end