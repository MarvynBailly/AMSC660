function [w, f, gnorm, k] = SG(fun,gfun, w, kmax, tol, bsz,n)
    f = zeros(kmax, 1);      
    gnorm = zeros(kmax, 1); 
    alpha = 0.3;
    %alpha = 0.01;

    for k = 1:kmax  
        % let's try some decreasing methods
        %alpha = 0.3/2^k;
        %alpha = 0.3/k^2;
        %alpha = 0.999*alpha;
        
        indices = randperm(n,bsz);
        g = gfun(indices, w);

        w = w - alpha * g;
        f(k) = fun(indices,w);
        gnorm(k) = norm(g);
        
        if gnorm(k) < tol
            break;
        end
        fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
    end
end
