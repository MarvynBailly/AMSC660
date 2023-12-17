    function question3()
    close all;

    N = 30;    
    iterMax = 1e8;
    betas = 0.2:0.01:1;
    
    meanMags = zeros(1,length(betas));
    varMags = zeros(1,length(betas));
    
    % Maybe vectorize the betas if time but for now for loop
    for i = 1:length(betas)
        % set up N x N mesh with all spins up
        mesh = ones(N,N);
        
        % use MCMC to find running mean and var
        [meanMag,varMag] = Metropolis(mesh,betas(i),iterMax);
        meanMags(i) = meanMag;
        varMags(i) = varMag;

        fprintf('Progress: %d%% complete: %d out of %d\n', round(i/length(betas)*100), i, length(betas));
    end

    %plot
    figure(1)
    grid on;


    plot(betas,meanMags,'ro', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    hold on;
    
    plot(betas, meanMags + sqrt(varMags),'--','Color', [1, 0.5, 0]);
    plot(betas, meanMags - sqrt(varMags),'--','Color', [1, 0.5, 0]);

    plotTrueMu();

    ylim([-0.5 1]);
    set(gca,'Fontsize',20);
    xlabel('beta','Fontsize',20);
    ylabel('Mean Magntization','Fontsize',20);

    figure(1)
    grid on;
    % save the figure
    saveas(figure(1), '..\images\coolplot3.png');
end

function plotTrueMu()
    betas =  linspace(0.2, 1, 1000);
    mu = zeros(size(betas));
    for i = 1:length(betas)
        if betas(i) > 0.4408
            mu(i) = (1 - sinh(2 * betas(i))^(-4)).^(1/8);
        else
            mu(i) = 0; 
        end
    end
    figure(1); 
    plot(betas, mu, 'b-', 'LineWidth', 2);
end

function [mu,var] = Metropolis(mesh,beta,iterMax)
    N = length(mesh);
    % compute magnetization
    m = sum(mesh,'all')/N^2;
    mu = m;
    var = 0;

    for iter = 1:iterMax
        % randomly pick a site, flip, and compute difference
        k = randi(N);
        l = randi(N);
        
        % compute delta H
        Dh = 2*mesh(k,l)*(mesh(mod(k-2, N) + 1, l) + mesh(mod(k, N) + 1, l) + mesh(k, mod(l-2, N) + 1) + mesh(k, mod(l, N) + 1));
        
        % if we accept this one, update the spins
        % note the rand is U(0,1)
        if Dh < 0 || rand < exp(-beta*Dh)
            mesh(k,l) = -mesh(k,l);
            % update the magnetization
            m = sum(mesh,'all')/N^2;
        end
    
        % update the running mean
        mu = (iter*mu+m)/(iter+1);
            
        % update the running variance
        var = ((iter - 1) * var + (m - mu)^2) / iter;
    end
end




