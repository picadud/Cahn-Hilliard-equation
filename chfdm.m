
function chfdm(N)

%default valuesch
space_length = 10;
dx = space_length / (N-1);
X = 0:dx:space_length;

%time_length = 100;
epsilon = 0.1;
M = 1;
dt = dx / 10;

% Initial condition
U = zeros(N,1);

%interior points
% for i = 3:N-2
%     U(i) = 0.1 * sin(2*pi*(i-3)*dx) ...
%           + 0.01 * cos(4*pi*(i-3)*dx) ...
%           + 0.06 * sin(4*pi*(i-3)*dx) ...
%           + 0.02 * cos(10*pi*(i-3)*dx);
% end
% 
% % boundary poidisnts
% U(1) = 0.03;
% U(2) = 0.03; 
% U(N-1) = 0.03; 
% U(N) = 0.03; 
%disp("U")
%disp(U)
for i = 1: N
    U(i) = 1.6 * i*dx ./  10  - 0.8;
end

%define 1d laplacian(central difference) matrix

% Create the main diagonal (with -2s)
main_diag = -2 * ones(N, 1);

% Create the diagonals above and below the main diagonal (with 1s)
upper_diag = ones(N-1, 1);
lower_diag = ones(N-1, 1);

% Assemble the matrix using the diag function
laplacian = (diag(main_diag) + diag(upper_diag, 1) + diag(lower_diag, -1));
%apply boundary condition to laplacian
laplacian(1,2) =  2;
laplacian(N,N-1) = 2;
%laplacian(1,N) = 1;
%laplacian(N, 1) = 1;
%display(laplacian)
squared_laplacian = laplacian^2;
%display(squared_laplacian)
laplacian = laplacian ./ dx^2;
%squared_laplacian = 


%free energy function
f = @(u) sum(epsilon.^2 ./ 2 .* gradient(u).^2 + (u.^2 -1).^2);

%composition plot
figure;

%plot(X, U, 'LineWidth', 2);

%free energy plot

%mass plot

t(1) = 0;
E(1) = f(U);
mass(1) = sum(U).\N;
i = 1;

%initiate explicit vector
U_exp = U;
while 1%i < 1000
    i = i+1;
    t(i) = (i-1) * dt;
    Un = diag(U);
    %assemble Eyre spliting matrix
    A = eye(N) + dt * (epsilon^2 * laplacian^2 - laplacian * Un^2);

    %assemble vector
    b =  U - dt .* laplacian * U;
    %slove the linear system
    U = A\b;
    %explicit solver
    U_exp = dt .* (-epsilon.^2 * laplacian^2*U_exp + laplacian * (U_exp.^3 - U_exp)) + U_exp;
    %disp(U_exp)
    %U = clip(U, -0.99, 0.99);
    %plot U
    subplot(2,2,1)
    plot(X, U, 'LineWidth', 2);
    ylim([-1,1])
    ylabel("composition u");
    xlabel("x");
    
    %plot U_exp
    subplot(2,2,2)
    plot(X, U_exp, 'LineWidth', 2);
    ylim([-1,1])
    ylabel("composition u_exp");
    xlabel("x");

    %plot free energy
    subplot(2,2,3)
    E(i) = f(U);
    plot(t, E)
    ylabel("total free energy");
    xlabel("t");
    pause(0.00001)
    
    %plot mass
    subplot(2,2,4)
    mass(i) = sum(U);
    plot(t, mass)
    xlabel("t");
    ylabel("mass");



end