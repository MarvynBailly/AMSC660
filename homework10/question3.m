function question3()
    
    rstar = 2^(1/6);
    Na = 7;
    model = 0;
    results = zeros(20,3);
    k = 1; 
    c = 0.25; 
    tol = 1e-6;

    for j = 1:10
        close all
        xyz = initial_configuration(model,Na,rstar);
       
        for i = 0:1
            %results(k,:) = LJ_line_search(i,c,xyz,Na)
            results(k,:) = LJ_trust_region(i,xyz);
            %LJ_trust_region(i,xyz);
            k = k + 1;
        end
    end
    results
    save('data.mat','results')
    figure(2)
    xlim([0,80])
    legend("Newton with exact","BFGS with Dogleg")
    figure(3)
    legend("Newton with exact","BFGS with Dogleg","Location","SouthEast")
end


%% make initial configuration
function xyz = initial_configuration(model,Na,rstar)
xyz = zeros(3,Na);
switch(model)
    case 1 % Pentagonal bipyramid
        p5 = 0.4*pi;
        he = sqrt(1 - (0.5/sin(0.5*p5))^2);
        for k = 1 : 5
            xyz(1,k) = cos((k-1)*p5); 
            xyz(2,k) = sin((k-1)*p5); 
            xyz(3,k) = 0;  
        end
        xyz(3,6) = he';
        xyz(3,7) = -he';

case 2 % Capped octahedron
        r = 1/sqrt(2);
        p4 = 0.5*pi;
        pp = p4/2;
        p0 = pi*1.5 - pp;
        x0 = sin(pp);
        y0 = cos(pp);
        z0 = 0;
        for k = 1 : 4
            xyz(1,k) = x0 + r*cos(p0 + (k-1)*p4);
            xyz(2,k) = y0 + r*sin(p0 + (k-1)*p4);
            xyz(3,k) = z0;
        end
        xyz(:,5) = [x0, y0, z0 + r]';
        xyz(:,6) = [x0, y0, z0 - r]';
        xyz(:,7) = [3*x0, y0, z0 + r]';

    case 3  % Tricapped tetrahedron
        p3 = 2*pi/3;
        pp = p3/2;
        r = 1/sqrt(3);
        beta = 0.5*pi -asin(1/3) - acos(1/sqrt(3));
        r1 = cos(beta);
        p0 = 1.5*pi - pp;
        x0 = sin(pp);
        y0 = cos(pp);
        z0 = 0;
        for k = 1 : 3
            xyz(1,k) = x0 + r*cos(p0 + (k-1)*p3);
            xyz(2,k) = y0 + r*sin(p0 + (k-1)*p3);
            xyz(3,k) = z0;
            xyz(1,k + 3) = x0 + r1*cos(p0 + pp + (k-1)*p3);
            xyz(2,k + 3) = y0 + r1*sin(p0 + pp + (k-1)*p3);
            xyz(3,k + 3) = z0 + sqrt(2/3) - sin(beta);
        end
        xyz(:,7) = [x0, y0, z0 + sqrt(2/3)]';

    case 4 % Bicapped trigonal bipyramid
        p3 = 2*pi/3;
        pp = p3/2;
        r = 1/sqrt(3);
        beta = 0.5*pi -asin(1/3) - acos(1/sqrt(3));
        r1 = cos(beta);
        p0 = 1.5*pi - pp;
        x0 = sin(pp);
        y0 = cos(pp);
        z0 = 0;
        for k = 1 : 3
            xyz(1,k) = x0 + r*cos(p0 + (k-1)*p3);
            xyz(2,k) = y0 + r*sin(p0 + (k-1)*p3);
            xyz(3,k) = z0;
        end
        xyz(:,4) = [x0 + r1*cos(p0 + pp), y0 + r1*sin(p0 + pp), z0 + sqrt(2/3) - sin(beta)]';
        xyz(:,5) = [x0 + r1*cos(p0 + pp + p3), y0 + r1*sin(p0 + pp+p3), z0 - sqrt(2/3) + sin(beta)]';
        xyz(:,6) = [x0, y0, z0 - sqrt(2/3)]';
        xyz(:,7) = [x0, y0, z0 + sqrt(2/3)]';

    otherwise % random configuration
        hR = 0.01;
        xyz = zeros(3,Na);
        xyz(:,1) = [0;0;0];
        a = randn(3,Na - 1);
        rad = sqrt(a(1,:).^2 + a(2,:).^2 + a(3,:).^2);
        a = a*diag(1./rad);
        for i = 2 : Na
            clear rad
            clear x
            rad = sqrt(xyz(1,1 : i - 1).^2 + xyz(2,1 : i - 1).^2 + xyz(3,1 : i - 1).^2);
            R = max(rad) + rstar;
            xa = R*a(:,i - 1);
            x = [xyz(:,1 : i - 1), xa];
            f = LJ(x(:));
            fnew = f;
            while 1
                R = R - hR;
                xa = R*a(:,i - 1);
                x = [xyz(:,1 : i - 1), xa];
                f = fnew;
                fnew = LJ(x(:));
                if fnew > f
                    break;
                end
            end
            xyz(:,i) = xa;
        end
        cmass = mean(xyz,2);
        xyz = xyz - cmass*ones(1,Na);
end
xyz = xyz*rstar;
end

function v = LJ(xyz) % 
global feval
feval = feval + 1;
m = length(xyz);
Na = m/3;
x = reshape(xyz,3,Na);
r = zeros(Na);
for k = 1 : Na
    r(k,:) = sqrt(sum((x - x(:,k)*ones(1,Na)).^2,1));
end
r = r + diag(ones(Na,1));
aux = 1./r.^6;
L = (aux - 1).*aux;
L = L - diag(diag(L));
v = 2*sum(sum(L));
end
