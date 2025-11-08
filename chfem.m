 

function chfem(N)
%dbstop if error  % Once an error occurs, execution stops; inspect workspace.
x = linspace(0,10,N);

dx = x(2)-x(1);
dt = dx ./ 100;
epsilon = 0.5;


% Initialise the matrix Ahat, Mhat,  and vector fhat.
A = sparse(N,N);
M = sparse(N,N);
f = zeros(N,1);

u_n = zeros(N,1);

% Initial condition
v_n = zeros(N,1);

%interior points
for i = 3:N-2
    v_n(i) = 0.1 * sin(2*pi*(i-3)*dx) ...
          + 0.01 * cos(4*pi*(i-3)*dx) ...
          + 0.06 * sin(4*pi*(i-3)*dx) ...
          + 0.02 * cos(10*pi*(i-3)*dx);
end
v_n(1) = 0.03;
v_n(2) = 0.03; 
v_n(N-1) = 0.03; 
v_n(N) = 0.03; 



% loop over the intervals
for k = 1:(length(x)-1)
    % extract endpoints of the interval
    x1 = x(k);
    x2 = x(k+1);
    % evaluate length of interval k.
    len = x2-x1;
    % evaluate derivatives of basisfunctions on interval k.
    dphi(1) = 1/(x1-x2);

    dphi(2) = 1/(x2-x1);

    %basis function and each end points
    basis = {
        @(y) 1 - (y-x1)./len
        @(y) (y-x1)./len
        };
    %concentration function for the current time n  defined on this interval and derivative of the double well potential c^3-c
    c = @(y) v_n(k) .* basis{1}(y) + v_n(k+1) .* basis{2}(y);
    f_c = @(y) c(y).^3 - c(y);
    
    % enumerate the basisfunctions on interval k.
    enum([1 2]) = [k k+1];
    for i=1:2
        % evaluate intergrals related to f.
        fmul1 = @(y) basis{i}(y) .* f_c(y);
        f(enum(i) ) = f(enum(i)) + integral(fmul1, x1,x2);
        for j=1:2
            % evaluate integral related to A, M
            A(enum(i),enum(j)) = A(enum(i),enum(j)) + dphi(i)*dphi(j)*len;
            fmul2 = @(y) basis{i}(y) .* basis{j}(y);
            M(enum(i),enum(j)) = M(enum(i),enum(j)) + integral(fmul2,x1,x2); 

        end
    end
end

% remove basisfunction 1 and N from the system



%while loop to calculate mu and c for t^n+1, assemble the global matrix and
%vector and solve
counter = 0;
figure
while counter <10000
    counter = counter + 1;
    plot(x, v_n)
    ylim([-1.5,1.5])
    pause(0.001)
    G = [M, dt.* A; -epsilon.^2 .* A, M];
    b = [M*v_n;  f];
    sol = G \ b; 
    v_n = sol(1:N);
    sprintf("mass: %f", sum(v_n)./N)
    u_n = sol(N+1:end);
    sprintf("energy: %f", sum(u_n))

    %reset the f vector
    f = zeros(N,1);
    for k = 1:(length(x)-1)
         % extract endpoints of the interval
        x1 = x(k);
        x2 = x(k+1);
        len = x2-x1;

        %basis function on each end point
        basis = {
        @(y) 1 - (y-x1)./len
        @(y) (y-x1)./len
        };
        c = @(y) v_n(k) .* basis{1}(y) + v_n(k+1) .* basis{2}(y);
        f_c = @(y) c(y).^3 - c(y);
        % enumerate the basisfunctions on interval k.
        enum([1 2]) = [k k+1];
        for i=1:2
            % evaluate intergrals related to f.
            fmul1 = @(y) basis{i}(y) .* f_c(y);
            f(enum(i) ) = f(enum(i)) + integral(fmul1, x1,x2);
        end
    end

end

     




% b = vhat(2:(N-1),1);
% u(1,1) = 0;
% u(N,1) = 0;
% u(2:(N-1),1) = A\b