% function [w, f, gnorm, k] = SG(fun,gfun, w, kmax, tol, bsz,n)
%     f = zeros(kmax, 1);      
%     gnorm = zeros(kmax, 1); 
%     alpha = 0.01;%3;
% 
%     for k = 1:kmax  % epochs
% 
%         indices = randperm(n,bsz);
%         g = gfun(indices, w);
%         alpha = 0.999*alpha;
%         w = w - alpha * g;
%         f(k) = fun(indices,w);
%         gnorm(k) = norm(g);
% 
% 
%         fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
%         % Check for convergence
%         if gnorm(k) < tol
%             break;
%         end
%     end
% end

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
        %fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
    end
    fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
end

function [w, f, gnorm, k] = SG111(fun,gfun, w, kmax, tol, bsz,n)
    f = zeros(kmax, 1);      
    gnorm = zeros(kmax, 1); 
    alpha = 0.3;

    for k = 1:kmax  % epochs
        randI = randperm(n);
        % let's try some decreasing methods
        %alpha = 0.3/2^k;
        %alpha = 0.3/k^2;
        alpha = 0.999*alpha;
        f_avg = zeros(ceil(n/bsz),1);
        gnorm_avg = zeros(ceil(n/bsz),1);
        %for i = 1:ceil(n/bsz)
        indices = randperm(n,bsz);
        %indices = randI((i-1)*bsz+1:min(i*bsz,n));
        g = gfun(indices, w);

        w = w - alpha * g;
        %f(k) = fun(indices,w);
        gnorm(k) = norm(g);
        %f_avg(i) = fun(indices,w);
        %gnorm_avg(i) = norm(g);
        % Check for convergence
        if gnorm(k) < tol
            break;
        end
        %end
        %divide k by bsz to get epoch
        %f(k) = fun(1:n,w);

        %f(k) = mean(f_avg);
        %g(k) = mean(gnorm_avg);

        fprintf('epoch = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
    end
end










function [w, f, gnorm, k] = SG110(fun,gfun, w, kmax, tol, bsz,n)
    f = zeros(kmax, 1);      
    gnorm = zeros(kmax, 1); 
    alpha = 0.3;

    for k = 1:kmax  % epochs
        randI = randperm(n);
        % let's try some decreasing methods
        %alpha = 0.3/2^k;
        %alpha = 0.3/k^2;
        alpha = 0.999*alpha;
        for i = 1:ceil(n/bsz)
            %randomBatch = randI(1:bsz);
            %indices = randperm(n,bsz);
            indices = randI((i-1)*bsz+1:min(i*bsz,n));

            %g = gfun(randomBatch,w);
            g = gfun(indices, w);
            w = w - alpha * g;
        end
        f(k) = fun(1:n,w);
        gnorm(k) = norm(gfun(1:n,w));

        fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
        % Check for convergence
        if gnorm(k) < tol
            break;
        end
    end
end


function [w, f, gnorm, k] = SG10(fun,gfun, w, kmax, tol, bsz,n)
    f = zeros(kmax, 1);      
    gnorm = zeros(kmax, 1); 
    alpha = 0.03;

    for k = 1:kmax  % epochs
        randI = randperm(n); 
        alpha = 0.999*alpha;
        for i = 1:ceil(n/bsz)
            %randomBatch = randI(1:bsz);
            indices = randperm(n,bsz);
            %indices = randI((i-1)*bsz+1:min(i*bsz,n));

            %g = gfun(randomBatch,w);
            g = gfun(indices, w);
            w = w - alpha * g;
        end
        f(k) = fun(n,w);
        gnorm(k) = norm(gfun(n,w));

        fprintf('k = %d, f = %d, alpha = %d, gnorm = %d\n',k,f(k),alpha, gnorm(k))
        % Check for convergence
        if gnorm(k) < tol
            break;
        end
    end
end













