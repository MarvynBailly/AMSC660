function question1()
    close all
    v = @(x)(exp(-(x^4 - 2*x^2 + 1)));

    a = 3;

    % Part a
    va = v(a);
    fprintf("v(%d) = %d\n",a,va);

    % check to see if a is big enough
    if va > 10e-16
        fprintf("Pick better a value\n");
        return;
    end

    Z = trapezoid(v,-a,a,500);
    fprintf("Z is approx %d\n",Z);

    % Part b
    [sigmaOpt,c] = findOptimalSigma(Z);
    fprintf("Found optimal sigma = %d and c = %d\n",sigmaOpt,c);

    % Part c
    f = @(x)(1/(Z) * exp(-(x.^4 - 2*x.^2 + 1)));
    g = @(x)(1/(sqrt(2*pi*sigmaOpt^2))*exp(-x.^2/(2*sigmaOpt^2)));
    samples = acceptReject(c,sigmaOpt,f,g);

    % Part d
    expectation = monteCarloInt(samples); 
    fprintf('E[|x|] = %f\n', expectation);

    % save figures 
    saveas(figure(1), 'images\q1max.png');
    saveas(figure(2), 'images\q1histo.png');
end

function expectation = monteCarloInt(samples)
    expectation = mean(abs(samples));
end

function eta = acceptReject(c,sigma,f,g)
    N = 1e8; % the number of samples
    v = randn(N,1);
    xi = sigma*v;
    u = rand(N,1);

    % compute f(xi)/(c g(xi))
    ind = find(u <= f(xi) ./ (c * g(xi)));
    Na = length(ind); % the number of accepted RVs
    eta = xi(ind);
    
    fprintf("N/Na = %d\n",N/Na);

    % plot a histogram to test the distribution
    nbins = 500; % the number of bins
    etamax = max(eta);
    etamin = min(eta);
    nb1 = nbins + 1;
    x = linspace(etamin,etamax,nb1);
    h = x(2) - x(1); % bin width
    % xc = centers of bins
    xc = linspace(etamin + 0.5*h,etamax - 0.5*h,nbins);
    hh = zeros(nbins,1); % heights of the bins
    for i = 1 : nbins
        ind = find(eta >= x(i) & eta < x(i + 1));
        hh(i) = length(ind);
    end
    hh = hh/(Na*h); % scale the histogram
    %f = g(x);
    f = f(x);
    figure(2);
    plot(x,f,'r','Linewidth',2);
    hold on;
    plot(xc,hh,'b','Linewidth',2);
    grid;
    set(gca,'Fontsize',20);
    xlabel('x','Fontsize',20);
    ylabel('f(x)','Fontsize',20);
    legend('True N(0,1)','Generated N(0,1)');
end

function [sigmaOpt,c] = findOptimalSigma(Z)
    syms x sigma
    
    % Function definitions
    f = 1/Z * exp(-(x^4 - 2*x^2 + 1));
    g = 1/sqrt(2*pi*sigma^2)*exp(-x^2/(2*sigma^2));
    
    % Ratio of f and g
    ratio = f / g;
    
    % Derivative of d with respect to x
    derRatio = diff(ratio, x);
    % Solve dd(x) == 0 for x and pick the third solution
    critPoints = solve(derRatio == 0, x);
    
    % see which one gives the max
    % first compute the second derivative
    secondDerRatio = diff(derRatio,x);
    
    % see which ones is negative
    figure(1)
    fplot(subs(secondDerRatio,x,critPoints),[-10,10])
    axis([-10 10 -200 50])
    legend("CP1","CP2","CP3")
    
    xStar = critPoints(3);
    
    % plug the max value back into ratio
    ratioMax = subs(ratio, x, xStar);
    
    % Then compute the min value 
    sigmaOpt = fminbnd(matlabFunction(ratioMax), -5, 5);
    c = double(subs(ratioMax,sigma,sigmaOpt));
end


function sol = trapezoid(fun,a,b,N)
    sol = 0;
    dx = (b-a)/N;
    approx = 0;
    for i = 0:N
        if i == 0 || i == N
            approx = fun(a);
        else
            approx = 2*fun(a + i*dx);
        end
        sol = sol + approx;
    end
    sol = dx/2 * sol;
end
