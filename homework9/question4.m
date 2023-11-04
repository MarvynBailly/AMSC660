function question4()
    close all;
    a = 100;
    func = @(x)(1-x(1)).^2 + a*(x(2) - x(1).^2).^2; 
    gfun = @(x)[-2*(1-x(1))-4*a*(x(2)-x(1)^2)*x(1);2*a*(x(2)-x(1)^2)]; 
    Hfun = @(x)[2 + 12*a*x(1)^2 - 4*a*x(2), -4*a*x(1); -4*a*x(1), 2*a];
    xstar = [1;1];
    %x0 = [1.2;1.2];
    x0 = [-1.2;1];
    
    for i = 1:5
        [alphaVals,norms,~,xvals,iter] = line_search(func,gfun,Hfun,xstar,x0,i,0);
        figure(2);
        %clf;
        hold on;
        grid on;
        semilogy(1:(iter-1),norms(2:iter),'Linewidth',2);
        xlabel('Iteration #');
        ylabel('||(xk,yk) - (x*,y*)||');
        ylim([10e-18 10e0])
        xlim([0 200])
        set(gca, 'YScale', 'log')
    end
    
    figure(2)
    linear_convergence = (0.5 .^ (1:(iter-1)));
    semilogy(1:(iter-1), linear_convergence, '--', 'LineWidth', 2);
    semilogy(96:(iter-1)+95, linear_convergence, '--', 'LineWidth', 2);
    legend("SD","N","BFGS","FRCG","PRCG","linear convergence","Location","southeast")

    % alpha values
    figure(1);
    %clf;
    %hold on;
    grid on;
    plot(0:(iter-1),alphaVals(1:iter),'Linewidth',2)
    xlabel('Iteration #');
    ylabel('Step size');
     
    % % norms
    % figure(2);
    % %clf;
    % hold on;
    % grid on;
    % semilogy(1:(iter-1),norms(2:iter),'Linewidth',2);
    % xlabel('Iteration #');
    % ylabel('||(xk,yk) - (x*,y*)||');
    % set(gca, 'YScale', 'log')
    % %linear_convergence = (0.5 .^ (1:(iter-1)));
    % %semilogy(1:(iter-1), linear_convergence, '--', 'LineWidth', 2);
    % %legend('Actual Convergence', 'Linear Convergence');

    % contour
    figure(3)
    [X,Y] = meshgrid(linspace(-2,2,400),linspace(-1,3,400));
    Z = (1-X).^2 + a*(Y - X.^2).^2;
    contour(X,Y,Z,logspace(-1,3,10), 'LineColor', 'k');
    xlim([-1.5 1.5]);
    ylim([-0.2 1.2]);
    xlabel('x');
    ylabel('y');
    hold on;
    plot(xvals(1:iter,1), xvals(1:iter,2), 'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    plot(xstar(1),xstar(2),'g.','Markersize',25,'MarkerFaceColor', 'g');

    
end

function [alphaVals,norms,fvals,xvals,iter] = line_search(func,gfun,Hfun,xstar,x0,dir,print_flag)
    tol = 1e-9; % stop iterations when || grad f|| < tol
    iter_max = 500; % the maximal number of iterations
    BFGSReset = 20;
    % parameters for backtracking line search
    c = 0.45;
    rho = 0.9;
    direction = dir;
    
    x = x0;
    f = func(x);
    g = gfun(x);
    normD = norm(x - xstar);
    %fprintf("Initially, f = %d, ||grad f|| = %d\n",f,norm_g);
    iter = 1;
    
    fvals = zeros(iter_max);
    fvals(1) = f;
    norms = zeros(iter_max);
    norms(1) = normD;
    alphaVals = zeros(iter_max);
    alphaVals(1) = 1;
    xvals = zeros(iter_max,2);
    xvals(1,:) = x;


    fail_flag = 0;
    
    B = eye(length(x));
    while normD > tol && iter < iter_max 
        % choose search direction
        switch direction
            case 1 % steepest descent
                p = -g;
                dir = "SD";
            case 2 % Newton
                H = Hfun(x);
                [~,flag] = chol(H);
                if flag == 0 % H is SPD, use Newton's direction
                    p = -H\g; 
                    dir = "Newton";
                else % use the steepest descent direction
                    p = -g;
                    dir = "SD";
                end
            case 3 % BFGS
                if iter == 1
                    p = -B\g;
                else
                    s = x - xtemp;
                    y = g - gtemp;
                    if(mod(iter,BFGSReset)==0)
                        B = eye(length(x),length(x));
                    else
                        B = B - (B*s*s'*B)/(s'*B*s) + (y*y')/(y'*s);
                    end
                    p = -B\g;
                end
                dir = "BFGS";
            case 4 % Fletcher-Reeves nonlinear CG
                if iter == 1
                    p = - g;
                else 
                    beta = (g'*g)/(gtemp'*gtemp);
                    p = -g + beta *p;
                end 
                dir = "FRCG";
            case 5 % Polak-Ribiere nonlinear CG
                if iter == 1
                    p = -g;
                else
                    beta = (g'*(g-gtemp))/norm(gtemp)^2;
                    beta = max(beta,0);
                    p = -g + beta*p;
                end
                dir = "PRCG";
            otherwise
                return
        end
        % normalize the search direction if its length greater than 1
        norm_p = norm(p);
        if norm_p > 1
            p = p/norm_p;
        end
        % do backtracking line search along the direction p
        a = 1;
        f_temp = func(x + a*p);
        cpg = c*p'*g;
        while f_temp > f + a*cpg % check Wolfe's condition 1
            a = a*rho;
            if a < 1e-14
                fprintf("line search failed\n");
                iter = iter_max;
                fail_flag = 1;
                break;
            end
            f_temp = func(x + a*p);        
        end
        xtemp = x;
        gtemp = g;
        
        x = x + a*p;
        f = func(x);
        g = gfun(x);
        
        normD = norm(x - xstar);
        
        if print_flag == 1
            fprintf("iter %d : dir = %s, f = %d, step length = %d, norm = %d\n",iter,dir,f,a,normD);
        end
        iter = iter + 1;
        fvals(iter) = f;
        norms(iter) = normD;
        alphaVals(iter) = a;
        
        xvals(iter,:) = x;
        if fail_flag == 1
            break;
        end
    end
    fprintf("iter %d : dir = %s, f = %d, step length = %d, norm = %d\n",iter,dir,f,a,normD);
end