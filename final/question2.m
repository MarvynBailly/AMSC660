function question2()
    close all;

    % let's define some things
    tol = 1e-4;
    iterMax = 1e5;

    % first we import the adjaceny matrix
    A = importdata('Adjacency_matrix.csv');
    
    N = length(A);

    % initial positions using Gaussian random variables with mean zero and standard deviation N
    xy = randn(2*N, 1);
    
    % define grad
    % notice that we are stacking x and y on top of each other
    grad = @(x,y)forces(x,y,A);

    % let's use adams to find the min
    [xy_min,gnorms] = adams(grad,xy,iterMax,tol,N);
    
    % the visual the dude
    plot_graph(xy_min(1:N),xy_min(N+1:end),A,1)

    % and let's make a nice plot of the norms
    figure(2)
    semilogy(1:length(gnorms),gnorms,'-o')
    set(gca,'Fontsize',20);
    xlabel('iteration','Fontsize',20);
    ylabel('gnorm','Fontsize',20);

    % save the figures
    saveas(figure(1), '..\images\coolplot1.png');
    saveas(figure(2), '..\images\coolplot2.png');
end



function [w,gnorm] = adams(grad, w, kmax, tol,N)
    gnorm = zeros(kmax, 1); 
    beta1 = 0.9;
    beta2 = 0.999; 
    e = 10e-8;
    alpha = 0.001;
    m = zeros(length(w),1);
    v = zeros(length(w),1);
    
    for k = 1:kmax
        g = -grad(w(1:N),w(N+1:end));
        gnorm(k) = norm(g);

        m = beta1 * m + (1 - beta1).*g;
        v = beta2 * v + (1 - beta2).*g.^2;

        mt = m ./ (1 - beta1^k);
        vt = v ./ (1 - beta2^k);
        
        w = w -  alpha*mt ./ (sqrt(vt)  + e);
        
        % Print progress every 5000 iteration
        if mod(k,5000)==0
            fprintf('k = %d, gnorm = %d\n',k,gnorm(k))
        end
        % Check for convergence
        if gnorm(k) < tol
            break;
        end
    end
    gnorm = gnorm(1:k);
    fprintf('Ended at k = %d with gnorm = %d\n',k,gnorm(k))
end


