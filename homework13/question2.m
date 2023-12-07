function question2()
    ds = [5, 10, 15, 20];
    %ds=[3];
    
    for i = 1:length(ds)
        fprintf("Setting d = %d\n",ds(i))
        % part a
        vol = estimateVolumeCubeIntersectionBall(ds(i),1e6);
        fprintf("a) vol = %d\n",vol)
        

        % part b
        vol = estimateVolumeBallIntersectionCube(ds(i),1e6);
        fprintf("b) vol = %d\n",vol)
    end
end

function vol = estimateVolumeBallIntersectionCube(d, n)
    % Generate n random points inside a d-dimensional unit ball
    xi = randn(n, d);
    radii = sqrt(sum(xi.^2, 2))*ones(1,d);
    % Scale points to be uniformly inside the ball
    scale = rand(n, 1).^(1/d) * ones(1,d);
    xi = (xi .* scale) ./ (radii);

    % Check if points are inside the unit cube
    ind  = find(max(abs(xi),[],2) <= 0.5);
    Na = length(ind);
    
    % Calculate the fraction of points inside the cube
    fractionInside = Na / n;

    % Volume of the d-dimensional unit ball
    ballVolume = pi^(d / 2) / (d / 2  * gamma(d / 2 ));
   
    % Approximate volume of the intersection
    vol = fractionInside * ballVolume;
end

function vol = estimateVolumeCubeIntersectionBall(d, numSamples)
    xi = rand(numSamples, d) - 0.5;
    
    ind = find(sum(xi.^2, 2) <= 1); %see if it's in the B^d
    Na = length(ind); 

    vol = Na / numSamples; % Volume of B^d \cap C^d
end