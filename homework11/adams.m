function [w, f, gnorm, k] = adams(fun,gfun, w, kmax, tol, bsz,n)
    f = zeros(kmax, 1);      
    gnorm = zeros(kmax, 1); 
    beta1 = 0.9;
    beta2 = 0.999; 
    e = 10e-8;
    alpha = 0.001;
    m = zeros(length(w),1);
    v = zeros(length(w),1);

    for k = 1:kmax 
        indices = randperm(n,bsz);
        
        g = gfun(indices, w);
        f(k) = fun(indices,w);
        gnorm(k) = norm(g);

        m = beta1 * m + (1 - beta1).*g;
        v = beta2 * v + (1 - beta2).*g.^2;

        mt = m ./ (1 - beta1^k);
        vt = v ./ (1 - beta2^k);
        
        w = w -  alpha*mt ./ (sqrt(vt)  + e);

        fprintf('k = %d, f = %d, gnorm = %d\n',k,f(k),gnorm(k))
        % Check for convergence
        if gnorm(k) < tol
            break;
        end
    end
end
