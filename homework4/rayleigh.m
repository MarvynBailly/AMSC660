function rayleigh(n,tol)
    A = rand(n);
    A = A' + A;
    %n = 100;
    v = rand(n,1);
    v = v/norm(v);
    k = 1;
    mu(k) = v'*A*v;
    
    %tol = 1e-12;
    I = eye(n);
    res = abs(norm(A*v - mu(k)*v)/mu(k));
    fprintf('k = %d: lam = %d\tres = %d\n',k,mu(k),res);
    while res > tol
        w = (A - mu(k)*I)\v;
        k = k + 1;
        v = w/norm(w);
        mu(k) = v'*A*v;
        res = abs(norm(A*v - mu(k)*v)/mu(k));
        fprintf('k = %d: lam = %d\tres = %d\n',k,mu(k),res);
    end
    
    errorArray = abs(mu(1:k-1) - mu(k)*ones(1,k-1));
    cubicErrorArray = errorArray(2:k-1)./errorArray(1:k-2).^3;
    x = linspace(1,k-2,k-2);
    plot(x,cubicErrorArray,'o-');
    xlabel("kth iteration")
    ylabel('$\frac{e_{k+1}}{e_k^3}$', 'interpreter', 'latex');
end