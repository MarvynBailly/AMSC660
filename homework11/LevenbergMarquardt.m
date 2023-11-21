function [w,f,gnorm,k] = LevenbergMarquardt(r_and_J,w,kmax,tol)
    Delta_max = 1;
    Delta_min = 1e-14;
    Delta = 1;

    rho_good = 0.75;
    rho_bad = 0.25;
    eta = 0.01;

    I = eye(length(w));
    
    [r,J] = r_and_J(w);
    g = J'*r;
    B = J'*J + (1e-6)*I; 
    norm_g = norm(g);

    f = zeros(kmax,1);
    gnorm = zeros(kmax,1);
    f(1) = 1/2 * norm(r)^2;
    gnorm(1) = norm_g;

    for k = 1:kmax
        if gnorm(k) < tol
            break;
        end

        pstar = -B\g; % unconstrained minimizer
        if norm(pstar) <= Delta
            p = pstar;
        else % solve constrained minimization problem
            lam = 1; % initial guess for lambda
            while 1 
                B1 = B + lam*I;
                C = chol(B1); % do Cholesky factorization of B
                p = -C\(C'\g); % solve B1*p = -g
                np = norm(p);
                dd = abs(np - Delta); % R is the trust region radius
    
                if dd < 1e-6
                    break
                end
    
                q = C'\p; % solve C^\top q = p
                nq = norm(q);
                lamnew = lam + (np/nq)^2*(np - Delta)/Delta;
    
                if lamnew < 0
                    lam = 0.5*lam;
                else
                    lam = lamnew;
                end
            end
        end

        % assess the progress
        wnew = w + p;
        [rnew,Jnew] = r_and_J(wnew);
        fnew = 1/2 * norm(rnew)^2;
        gnew = Jnew'*rnew;
        mnew = f(k) + g'*p + 0.5*p'*B*p;
        rho = (f(k) - fnew+1e-14)/(f(k) - mnew+1e-14);
        % adjust the trust region
        if rho < rho_bad
            Delta = max([0.25*Delta,Delta_min]);
        else
            if rho > rho_good && norm(p) == Delta
                Delta = min([Delta_max,2*Delta]);
            end
        end
    
        % accept or reject step
        if rho > eta  
            w = wnew;
            g = gnew;
            r = rnew; 
            B = Jnew'*Jnew + (1e-6)*I;
            %fprintf("Accepted: iter # %d: f = %.10f, |df| = %.4e,rho = %.4e, Delta = %.4e\n",k,f,norm_g,rho,Delta);
            %fprintf('Accept: iter # %d: f = %.10f, |df| = %.4e, rho = %.4e, Delta = %.4e, j_sub = %d\n',iter,r,norm_g,rho,Delta,j_sub);
        else
            %fprintf("Rejected: iter # %d",iter);
            %fprintf("Accepted: iter # %d: f = %.10f, |df| = %.4e,rho = %.4e, Delta = %.4e\n",k,f,norm_g,rho,Delta);
            %fprintf('Reject: iter # %d: f = %.10f, |df| = %.4e, rho = %.4e, Delta = %.4e, j_sub = %d\n',iter,r,norm_g,rho,Delta,j_sub);
        end
        f(k+1) = 1/2 * norm(r)^2;
        gnorm(k+1) = norm(g);
    end
end