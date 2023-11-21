function [w, f, gnorm, k] = nesterov(fun,gfun, w, kmax, tol, bsz,n)
    f = zeros(kmax, 1);      
    gnorm = zeros(kmax, 1); 
    alpha = 0.01;
    

    y = w;
    wtemp = w;

    for k = 1:kmax  % epochs
        indices = randperm(n,bsz);
        mu = 1 - 3 / (5 + k);
            
        y = (1 + mu)*w - mu*wtemp;
        
        f(k) = fun(indices,y);
        g = gfun(indices, y);
        
        wtemp = w;
        w = y - alpha*g;
        g = gfun(indices, w);
        gnorm(k) = norm(g);


        %fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
        % Check for convergence
        if gnorm(k) < tol
            return;
        end
    end
    fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
end


function [w, f, gnorm, k] = nesterov1(fun,gfun, w, kmax, tol, bsz,n)
    f = zeros(kmax, 1);      
    gnorm = zeros(kmax, 1); 
    alpha = 0.3;
    

    y = w;
    wtemp = w;

    for k = 1:kmax  % epochs
        randI = randperm(n);
        for i = 1:ceil(n/bsz)
            %indices = randperm(n,bsz);
            indices = randI((i-1)*bsz+1:min(i*bsz,n));

            % y = w - alpha * g;
            % lambdaOld = lambda;
            % lambda = 1/2 * (sqrt(1 + 4*lambda^2));
            % gamma = (1 - lambdaOld) / (lambda);
            % w = (1 - gamma)*y + gamma*y;
            
            mu = 1 - 3 / (5 + k);
            
            y = (1 + mu)*w - mu*wtemp;
            
            f(k) = fun(indices,y);
            g = gfun(indices, y);
            
            wtemp = w;
            w = y - alpha*g;
            g = gfun(indices, w);
            gnorm(k) = norm(g);
        end


        fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
        % Check for convergence
        if gnorm(k) < tol
            break;
        end
    end
end