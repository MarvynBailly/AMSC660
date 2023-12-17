function f = forces(x,y,A)
% f = force = - grad U = column vector with 2*N components
% x, y are column vectors with N components
% A is an N-by-N adjacency matrix

N = length(x);
%% find pairwise distances between linked vertices
xaux = x*ones(size(x))';
yaux = y*ones(size(y))';
dx = A.*xaux - A.*(xaux');
dy = A.*yaux - A.*(yaux');
dxy = sqrt(dx.^2+dy.^2);

%% spring forces due to linked vertices
Aind = find(A == 1);
idiff = zeros(N);
idiff(Aind) = 1 - 1./dxy(Aind);
fx = -sum(idiff.*dx,2);
afx = min(abs(fx),1);
sfx = sign(fx);
fx = afx.*sfx;

fy = -sum(idiff.*dy,2);
afy = min(abs(fy),1);
sfy = sign(fy);
fy = afy.*sfy;

f_linked = [fx;fy];

% Calculate energy for linked vertices
%U_linked = 0.5 * sum((dxy(Aind) - 1).^2);

%% repelling spring forces due to unlinked vertices
h = sqrt(3);
Aind = find(A==0);
A = ones(size(A))-A;
dx = A.*xaux - A.*(xaux'); 
dy = A.*yaux - A.*(yaux');
dxy = sqrt(dx.^2+dy.^2);
fac = zeros(N);
diff = dxy - h;
fac(Aind) = min(diff(Aind),0); 
fx = sum(fac.*dx,2);
fy = sum(fac.*dy,2);
f_unlinked = -[fx;fy];

% Calculate energy for unlinked vertices
%U_unlinked = 0.5 * sum(min(dxy(Aind) - h, 0).^2);

f = f_linked + f_unlinked;
%U = U_linked + U_unlinked;
end
