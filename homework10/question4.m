function question4()
    close all;
    a = 100;
    func = @(x)(1-x(1)).^2 + a*(x(2) - x(1).^2).^2; 
    gfun = @(x)[-2*(1-x(1))-4*a*(x(2)-x(1)^2)*x(1);2*a*(x(2)-x(1)^2)]; 
    Hfun = @(x)[2 + 12*a*x(1)^2 - 4*a*x(2), -4*a*x(1); -4*a*x(1), 2*a];
    xstar = [1;1];
    %x0 = [1.2;1.2];
    x0 = [-1.2;1];
    
    for i = 0:1
        [norms,~,xvals,iter] = trust_region(func,gfun,Hfun,xstar,x0,i);
        figure(3);
        %clf;
        hold on;
        grid on;
        semilogy(1:(iter-1),norms(2:iter),'Linewidth',2);
        xlabel('Iteration #');
        ylabel('||(xk,yk) - (x*,y*)||');
        ylim([10e-10 50])
        xlim([0 400])
        set(gca, 'YScale', 'log')
        legend("Newton with exact","BFGS with Dogleg");

        % contour
        figure()
        [X,Y] = meshgrid(linspace(-2,2,400),linspace(-1,3,400));
        Z = (1-X).^2 + a*(Y - X.^2).^2;
        contour(X,Y,Z,logspace(-1,3,10), 'LineColor', 'k');
        %xlim([0.9,1.3]);
        %ylim([0.9,1.5]);
        xlim([-1.5 1.5]);
        ylim([-0.2 1.2]);
        xlabel('x');
        ylabel('y');
        hold on;
        plot(xvals(1:iter,1), xvals(1:iter,2), 'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        plot(xstar(1),xstar(2),'g.','Markersize',25,'MarkerFaceColor', 'g');
    end
    figure(1) 
    fig = gca;
    name = "4-cont-N-"+strjoin(string(x0), ',')+".png";
    exportgraphics(fig,"images/"+name)
    figure(2) 
    fig = gca;
    name = "4-cont-BFGS-"+strjoin(string(x0), ',')+".png";
    exportgraphics(fig,"images/"+name)
    figure(3) 
    fig = gca;
    name = "4-norm"+strjoin(string(x0), ',')+".png";
    exportgraphics(fig,"images/"+name)
end

function [norms,fvals,xvals,iter] = trust_region(func,gfun,Hfun,xstar,x0,direction)
    tol = 1e-9; % stop iterations when || grad f|| < tol
    iter_max = 400; % the maximal number of iterations
    BFGSReset = 20;

    % parameters for trust region
    Delta_max = 5; % the max trust-region radius
    Delta_min = 1e-12; % the minimal trust-region radius
    Delta = 1; % the initial radius
    eta = 0.1; % step rejection parameter
    subproblem_iter_max = 5; % the max # of iteration for quadratic subproblems
    tol_sub = 1e-1; % relative tolerance for the subproblem
    rho_good = 0.75;
    rho_bad = 0.25;

    x = x0;
    f = func(x);
    g = gfun(x);
    normD = norm(x - xstar);
    fprintf("Initially, f = %d, ||grad f|| = %d\n",f,normD);
    iter = 1;
    
    fvals = zeros(iter_max);
    fvals(1) = f;
    norms = zeros(iter_max);
    norms(1) = normD;
    xvals = zeros(iter_max,2);
    xvals(1,:) = x;

    I = eye(length(x));
    B = eye(length(x));
    
    while normD > tol && iter < iter_max   
        % solve the constrained minimization problem 
        if (direction == 0)
            B = Hfun(x);
        end
    
        flag_boundary = 0;
        % check if B is SPD
        eval_min = min(eig(B));
        j_sub = 0;
        if eval_min > 0 % B is SPD: B = R'*R, R'*R*p = -g 
            if direction == 0
                p = -B\g;
                p_norm = norm(p);
                if p_norm > Delta % else: we are done with solbing the subproblem
                    flag_boundary = 1;
                end
            elseif direction == 1 % use Dog leg method
                pu = - (norm(g)^2*g)/(g'*B*g);
                pu_norm = norm(pu);
                if pu_norm >= Delta
                    fprintf("pu >= delta\n");
                    tau = Delta/(pu_norm);
                    p = tau * pu;
                    flag_boundary = 1;
                else 
                    pb = -B\g;
                    if pb <= Delta
                        fprintf("pb <= delta\n");
                        p = pb;
                    else
                        % solve the quadratic in alpha
                        a = (pu - pb)'*(pu - pb);
                        b = 2*(pu'*(pb - pu));
                        c = pu'*pu - Delta^2;
                        
                        % get alpha
                        fprintf("pu < delta\n");
                        alpha = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
                        p = pu + alpha*(pb - pu);
                        flag_boundary = 1;
                    end
                end
            end
        else
            flag_boundary = 1;
        end
        if flag_boundary == 1 % solution lies on the boundary
            lambda_min = max(-eval_min,0);
            lambda = lambda_min + 1;
            R = chol(B+lambda*I);
            flag_subproblem_success = 0;
            while j_sub < subproblem_iter_max 
                j_sub = j_sub + 1; 
                p = -R\(R'\g);
                p_norm = norm(p);
                dd = abs(p_norm - Delta);
                if dd < tol_sub*Delta
                    flag_subproblem_success = 1;
                    break
                end
                q = R'\p;
                q_norm = norm(q);
                dlambda = ((p_norm/q_norm)^2)*(p_norm - Delta)/Delta;
                lambda_new = lambda + dlambda;
                if lambda_new > lambda_min
                    lambda = lambda_new;
                else
                    lambda = 0.5*(lambda + lambda_min);
                end
                R = chol(B+lambda*I);
            end
            if flag_subproblem_success == 0
                p = cauchy_point(B,g,Delta);
            end
        end
        % assess the progress
        
        xnew = x + p;
        fnew = func(xnew);
        gnew = gfun(xnew);
        mnew = f + g'*p + 0.5*p'*B*p;
        rho = (f - fnew+1e-14)/(f - mnew+1e-14);
        % adjust the trust region
        if rho < rho_bad
            Delta = max([0.25*Delta,Delta_min]);
        else
            if rho > rho_good && flag_boundary == 1
                Delta = min([Delta_max,2*Delta]);
            end
        end
        
        if direction == 1 && rho > eta % update bfgs matrix B if accepted
            p = xnew - x;
            y = gnew - g;
            if mod(iter,BFGSReset)==0
                B = eye(length(x));
            else
                B = B-((B*p)*(p'*B'))/(p'*B*p) +(y*y')/(y'*p);
            end
        end
    
        % accept or reject step
        if rho > eta  
            x = xnew;
            f = fnew;
            g = gnew;
            normD = norm(g);
            fprintf('Accept: iter # %d: f = %.10f, |(x,y) - (x^*,y^*)| = %.4e, rho = %.4e, Delta = %.4e, j_sub = %d\n',iter,f,normD,rho,Delta,j_sub);
        else
            fprintf('Reject: iter # %d: f = %.10f, |(x,y) - (x^*,y^*)| = %.4e, rho = %.4e, Delta = %.4e, j_sub = %d\n',iter,f,normD,rho,Delta,j_sub);
        end
        iter = iter + 1;
        fvals(iter) = f;
        norms(iter) = normD;
        xvals(iter,:) = x;
    end
end