function question2()
    data = readtable('MovieRankingData.csv','NumHeaderLines',1);
    % to get submatrix, take first eight columns and delete rows with
    % missing values
    dataMatrix = data{:,2:9};
    dataMatrix(any(isnan(dataMatrix), 2), :) = [];

    tol = 10e-12;
    k = 7;
        
    % by visual inspection
    completeSubMatrix = dataMatrix(1:2,:);

    % part a - project 
    errsPJD = PJD(completeSubMatrix,k,tol);
    
    % part b - Lee-Seung
    errsLS = leeseung(completeSubMatrix,k,tol);

    % plot
    figure(1);
    hold on;
    plot(0:length(errsPJD)-1,errsPJD);
    plot(0:length(errsLS)-1,errsLS);
    hold off;
    legend({"Forbenius norm squared of PJD","Forbenius norm squared of LS"});
    xlabel('Number of Iterations');
    ylabel('Forbineus norm squared');
    axis([0, length(errsPJD), 0, max(max(errsPJD),max(errsLS))])
    % plot
    figure(2);
    hold on;
    plot(0:length(errsPJD)-1,errsPJD);
    plot(0:length(errsLS)-1,errsLS);
    hold off;
    legend({"Forbenius norm squared of PJD","Forbenius norm squared of LS"});
    xlabel('Number of Iterations');
    ylabel('Forbineus norm squared');
    axis([0, length(errsLS), 0, max(max(errsPJD),max(errsLS))])
end

function errs = PJD(A,k,tol)
    % input: complete matrix A nxd
    % output: W nxk and H kxd
    
    % compute stepsize
    a = 1/100;

    % compute dims
    [n,d] = size(A);
    
    % define W,H,R
    W = rand(n,k);
    H = rand(k,d);
    R = A - W*H;
    errs = [norm(R,'fro')]
    % step away
    while tol < norm(R)
        W = max(W + a*R*H',0);
        R = A - W*H;
        H = max(H + a*W'*R,0);
        R = A - W*H;
        errs = [errs, norm(R,'fro')];
    end
end

function errs = leeseung(A,k,tol)
    % compute dims
    [n,d] = size(A);
    
    % define W,H,R
    W = rand(n,k);
    H = rand(k,d);
    R = A - W*H;
    errs = [norm(R,'fro')]

    while tol < norm(R)
        W = (W .* (A*H')) ./ (W*(H*H'));
        H = (H .* (W'*A)) ./ ((W'*W)*H);
        R = A - W*H;
        errs = [errs, norm(R,'fro')];
    end
end