function NestedDissection_anysize()
close all
fsz = 20; % fontsize
n = 20; % the set of unknown grid points is n-by-n
%[A,b,sol] = TestMatrixA(n + 2);

%A,b,sol from maze question:
load('A.mat');
load('b.mat');
%load('sol.mat');




A = -A;
b = -b;

%fix the A matrix
Atemp =  zeros(length(A)+2,length(A)+2);
Atemp(2:end-1,2:end-1) = A;
Atemp(1,1) = 1;
Atemp(end,end) = 1;
A = Atemp;

%fix the b vector
btemp = zeros(length(b)+2,1);
btemp(2:end-1) = b;
btemp(1) = 0;
btemp(end) = 1;
b = btemp; 

sol = A\b;


figure(1000)
spy(A)

%% Nested Dissection
% Nested Dissection algorithm solves Ax=b
% It rearranges A using a permutation matrix P that it constructs
% Note: P' = P^{-1}
% Hence it solves PAP'Px = Pb
% It decomposes PAP' = LU
% Hence we need to solve LUPx = Pb
% It first solves Ly = Pb by y = L\(Pb)
% Then it solves UPx = y by x = P'(U\y)

%%
% The grid size is n-by-n
level = 0;
[L,U,P,A] = MyDissection(A,n,n,level);
y = L\(P*b);
x = P'*(U\y);
fprintf('norm(x - sol) = %d\n',norm(x - sol));
figure;
spy(A);
set(gca,'fontsize',fsz);
grid
title(sprintf('Sparsity pattern of P*A*P^T for nx = %d, ny = %d',n,n),'Fontsize',20);
axis ij

end

%%
function [L,U,P,A] = MyDissection(A,nx_grid,ny_grid,level)
A0 = A;
[m,n] = size(A);
if m ~= n
    fprintf("A is not square: (%d-by-%d)\n",m,n);
    L=0;
    U=0;
    P=0;
    A=0;
    return
end
% if level is even do vertical split
% if level is odd do horizontal split
N_grid = nx_grid*ny_grid; % the grid size              % one flop here
par = mod(level,2); % parity                           % one flop here
switch par
    case 0 % vertical split
        if nx_grid >= 3
            nx_Omega1 = floor(nx_grid/2); % # of columns in Omega1  
            N_Omega1 = nx_Omega1*ny_grid; % |Omega1| 
            ind_Omega3 = N_Omega1 + 1 : N_Omega1 + ny_grid; % indices of Omega3
            N_Omega3 = length(ind_Omega3);
            ind_Omega1 = 1 : N_Omega1; % indices of Omega1
            ind_Omega2 = N_Omega1 + N_Omega3 + 1 : N_grid; % indices of Omega2
            N_Omega2 = length(ind_Omega2);
            ny_Omega1 = ny_grid;
            nx_Omega2 = nx_grid-nx_Omega1-1;
            ny_Omega2 = ny_grid; 
        else
            %[L,U] = lu(A);
            R = chol(A);
            U = R;
            L = R';

            P = speye(N_grid);
            return
        end    
    case 1 % horizontal split
        if ny_grid >= 3
            ny_Omega1 = floor(ny_grid/2); % # of rows in Omega1
            N_Omega1 = ny_Omega1*nx_grid; % |Omega1|
            ind_Omega3 = ny_Omega1 + 1 : ny_grid : N_grid; % indices of Omega3
            N_Omega3 = length(ind_Omega3);
            [ii,jj] = meshgrid(1 : ny_Omega1,1 : nx_grid);
            ind_Omega1 = sort(sub2ind([ny_grid,nx_grid],ii(:)',jj(:)'),'ascend'); % indices of Omega1
            [ii,jj] = meshgrid(ny_Omega1 + 2 : ny_grid,1 : nx_grid);
            ind_Omega2 = sort(sub2ind([ny_grid,nx_grid],ii(:)',jj(:)'),'ascend'); % indices of Omega2 
            N_Omega2 = length(ind_Omega2);
            nx_Omega1 = nx_grid;
            nx_Omega2 = nx_grid;
            ny_Omega2 = ny_grid-ny_Omega1-1; 
        else
            %[L,U] = lu(A);
            R = chol(A);
            U = R;
            L = R';

            P = speye(N_grid);
            return
        end    
    otherwise
        fprintf('Error: par = %d\n',par);
        return
end
%fprintf('size(A) = [%d,%d], nxy = %d, N_Omega1 = %d, N_Omega2 = %d, N_Omega3 = %d\n',size(A,1),size(A,2),ny_grid,N_grid,N_Omega1,N_Omega2,N_Omega3);
%fprintf('ind_Omega2(1) = %d, ind_Omega2(end) = %d\n',ind_Omega2(1),ind_Omega2(end));

A11 = A(ind_Omega1,ind_Omega1);
A22 = A(ind_Omega2,ind_Omega2);
[L11,U11,P11,~] = MyDissection(A11,nx_Omega1,ny_Omega1,level + 1);
[L22,U22,P22,~] = MyDissection(A22,nx_Omega2,ny_Omega2,level + 1);

P1 = speye(N_grid);
P1(ind_Omega1,ind_Omega1) = P11;
P1(ind_Omega2,ind_Omega2) = P22;
% set up the permutation matrix P
% this command puts ones in positions 
% with indices (1 : nxy,[ind1(:)',ind2(:)',ind(:)'])
P = sparse(1 : N_grid,[ind_Omega1(:)',ind_Omega2(:)',ind_Omega3(:)'],ones(1,N_grid));
P = P*P1;
A = P*A0*P';
% extract nonzero blocks of A
A11 = A(1 : N_Omega1,1 : N_Omega1);

istart2 = N_Omega1 + 1;
ifinish2 = N_Omega1 + N_Omega2;
istart3 = ifinish2 + 1;

A22 = A(istart2 : ifinish2,istart2 : ifinish2);
A13 = A(1 : N_Omega1,istart3 : end);
A23 = A(istart2 : ifinish2,istart3 : end);
A31 = A(istart3 : end,1 : N_Omega1);
A32 = A(istart3 : end,istart2 : ifinish2);
A33 = A(istart3 : end,istart3 : end);
% compute the Schur complement
%disp("----Computing the Schur Complement----")
fprintf('Size of A31 %dx%d Size of U11 and L11 %dx%d.\n',size(A33,1),size(A33,2),size(U11,1),size(U11,2))
fprintf('Size of A32 %dx%d Size of U22 and L22 %dx%d.\n',size(A32,1),size(A32,2),size(U22,1),size(U22,2))

figure(3)
spy(L11)
title(sprintf('Sparsity pattern of L11'));
figure(4)
spy(U11)
title(sprintf('Sparsity pattern of U11'));
figure(5)
spy(A13)
title(sprintf('Sparsity pattern of A13'));
figure(6)
spy(A31)
title(sprintf('Sparsity pattern of A31'));
figure(7)
spy(L22)
title(sprintf('Sparsity pattern of L22'));
figure(8)
spy(U22)
title(sprintf('Sparsity pattern of U22'));
figure(9)
spy(A23)
title(sprintf('Sparsity pattern of A23'));
figure(10)
spy(A32)
title(sprintf('Sparsity pattern of A32'));

S33 = A33 - A31*(U11\(L11\A13)) - A32*(U22\(L22\A23));



% compute LU factorization of S33
%[L33,U33] =     (S33);
R33 = chol(S33);
U33 = R33;
L33 = R33';

% form the LU decomposition of A
L = sparse(N_grid,N_grid);
L(1 : N_Omega1,1 : N_Omega1) = L11;
L(istart2 : ifinish2,istart2 : ifinish2) = L22;
L(istart3 : end,istart3 : end) = L33;
L(istart3 : end,1 : N_Omega1) = (U11'\A31')';
L(istart3 : end,istart2 : ifinish2) = (U22'\A32')';
U = sparse(N_grid,N_grid);
U(1 : N_Omega1,1 : N_Omega1) = U11;
U(istart2 : ifinish2,istart2 : ifinish2) = U22;
U(istart3 : end,istart3 : end) = U33;
U(1 : N_Omega1,istart3 : end) = L11\A13;
U(istart2 : ifinish2,istart3 : end) = L22\A23;

% fprintf('nx = %d, ny = %d, level = %d, size(A) = [%d,%d]\n',nx_grid,ny_grid,level,n,n);
% fprintf('norm(full(A - L*U)) = %d\n',norm(full(A - L*U)));
end



        
        
        
        
        

        
