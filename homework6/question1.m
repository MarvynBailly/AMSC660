function question1()
    data = readtable('MovieRankingData.csv','NumHeaderLines',1);
    dataMatrix = data{:,2:end};
    % save NaN locations
    nanLocations = isnan(dataMatrix);
    % make NaN zeros
    dataMatrix(nanLocations) = 0;
    % define things
    lambdas = [0.01,0.1,1,10,100];
    ks = 1:7;

    fprintf('My row before completion:[%s]\n', join(string(dataMatrix(79,:)), ','));

    % part a - 
    fprintf('My row after Low-Rank factorization:\n')
    for i = 1:length(lambdas)
        for j = 1:length(ks)
            [X,Y] = lowRankFact(dataMatrix,lambdas(i),ks(j),nanLocations);
            A = X*Y;
            fprintf('lambda = %d and k = %d:[%s]\n',lambdas(i),ks(j),join(string(A(79,:)), ','));
            fprintf('lambda = %d and k = %d:[%s]\n',lambdas(i),ks(j),join(string(A(79,:) - dataMatrix(79,1:10)), ','));
            %fprintf('lambda = %d and k = %d:[%s]\n',lambdas(i),ks(j),join(string(round(abs(A(79,1:10)),5)), ','));
            %fprintf('lambda = %d and k = %d:[%s]\n',lambdas(i),ks(j),join(string(round(abs(A(79,1:10) - dataMatrix(79,1:10)),1)), ','));
        end
    end
    % part b -
    for i = 1:length(lambdas)
        A = nuclear(dataMatrix,lambdas(i),nanLocations);
        fprintf('lambda = %d:[%s]\n',lambdas(i),join(string(round(abs(A(79,1:10)),1)), ','));
    end

    % test times
    tic
    [X,Y] = lowRankFact(dataMatrix,lambdas(1),ks(1),nanLocations);
    toc
    tic
    A = nuclear(dataMatrix,lambdas(1),nanLocations);
    toc
end

function [X,Yt] = lowRankFact(A, lambda,k,nanLocations)
    % input: incomplete matrix A
    % output: complete matrix X,Y 
    
    % get size of A
    [n,d] = size(A);

    % define X and Y^T
    X = rand(n,k);
    Yt = rand(k,d);

    for m = 1:1000
        for i = 1:n
            % compute Y_{omega_i}
            YOmega = Yt;
            YOmega(:,nanLocations(i,:)) = [];
            % compute a_{omega_i}
            bt = A(i,:);
            bt(nanLocations(i,:)) = [];
            % update
            HHt = YOmega*YOmega';
            X(i,:) = (HHt + lambda*eye(size(HHt))) \ YOmega*bt';
        end
        for j = 1:d
            % compute X_{omega_j} 
            XOmega = X;
            XOmega(nanLocations(:,j),:) = [];
            % compute a_{omega_i}
            bt = A(:,j);
            bt(nanLocations(:,j)) = [];
            % update 
            WtW = XOmega'*XOmega;
            Yt(:,j) = (WtW+lambda*eye(size(WtW))) \ XOmega'*bt;
        end
    end
    
end

function M = nuclear(A,lambda,nanLocations)
    %input: incomplete matrix A
    %output: complete matrix B
    M = rand(size(A));

   for i=1:1000
        % get A - M
        AM = A - M;
        % apply P_Omega
        AM(nanLocations) = 0;
        Pomega = AM;
        % get SVD of A - use econ for faster computation times
        [U,S,V] = svd(M + Pomega,'econ');
        % compute sLambda 
        Slambda = max(S - lambda*eye(size(S)), 0);
        % compute SLambda(A) and update M
        M = U*Slambda*V';
    end
end
