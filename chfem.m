 

function chfem(N, splitting)
%dbstop if error  % Once an error occurs, execution stops; inspect workspace.
%this script has quite bad style as for duplicate code, but since it is
%only for study purpose so please forgive me for that.
x = linspace(0,10,N);

dx = x(2)-x(1);
dt = dx ./ 10;
epsilon = 0.1;


% Initialise the matrix Ahat, Mhat,  and vector fhat.
A = sparse(N,N);
M = sparse(N,N);
F_c = sparse(N,N); %linearized contractive term
f_c_c = zeros(N,1);
f = zeros(N, 1);

u_n = zeros(N,1);

% Initial condition
v_n = zeros(N,1);

%interior points
% for i = 3:N-2
%     v_n(i) = 0.1 * sin(2*pi*(i-3)*dx) ...
%           + 0.01 * cos(4*pi*(i-3)*dx) ...
%           + 0.06 * sin(4*pi*(i-3)*dx) ...
%           + 0.02 * cos(10*pi*(i-3)*dx);
% end
% v_n(1) = 0.03;
% v_n(2) = 0.03; 
% v_n(N-1) = 0.03; 
% v_n(N) = 0.03; 
for i = 1: N
    v_n(i) = 1.6 * i*dx ./  10  - 0.8;
end
    


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
    f_c_e = @(y) - c(y); %for convex splitting
    f_all = @(y) c(y).^3 - c(y); % for other splitting

    
    % enumerate the basisfunctions on interval k.
    enum([1 2]) = [k k+1];
    for i=1:2
        % evaluate intergrals related to f.
        fmul1 = @(y) basis{i}(y) .* f_c_e(y);
        f_c_c(enum(i) ) = f_c_c(enum(i)) + integral(fmul1, x1,x2);
        fmul11 = @(y) basis{i}(y) .* f_all(y);
        f(enum(i)) = f_c_c(enum(i)) + integral(fmul11, x1,x2); 
        for j=1:2
            % evaluate integral related to A, M
            A(enum(i),enum(j)) = A(enum(i),enum(j)) + dphi(i)*dphi(j)*len;
            fmul2 = @(y) basis{i}(y) .* basis{j}(y);
            M(enum(i),enum(j)) = M(enum(i),enum(j)) + integral(fmul2,x1,x2); 
            fmul3 = @(y) c(y).^2 .* basis{i}(y) .* basis{j}(y); 
            F_c(enum(i),enum(j)) = F_c(enum(i),enum(j)) + integral(fmul3,x1,x2);
            

        end
    end
end

% remove basisfunction 1 and N from the system

%free energy function
free_engergy = @(u) sum(epsilon.^2 ./ 2 .* gradient(u).^2 + (u.^2 -1).^2);

%while loop to calculate mu and c for t^n+1, assemble the global matrix and
%vector and solve
t(1) = 0;
E(1) = 0;
mass(1) = sum(v_n)./N;
counter = 1;
figure
if splitting == 1
    while counter <10000
        counter = counter + 1;

        subplot(2,2,1)
        plot(x, v_n)
        ylim([-1.5,1.5])
        pause(0.00001)

        G = [M, dt.* A; -epsilon.^2.* A - F_c, M];
        b = [M*v_n;  f_c_c];
        sol = G \ b; 
        v_n = sol(1:N);
        ylabel("composition u");
        xlabel("x");

        t(counter) = t(counter-1) + dt;

        subplot(2,2,2)
        mass(counter) = sum(v_n)./N;
        plot(t, mass)
        u_n = sol(N+1:end);
        xlabel("t");
        ylabel("mass");
        

        subplot(2,2,3)
        E(counter) = free_engergy(v_n);
        plot(t, E)
        ylabel("total free energy");
        xlabel("t");
        
    
        %reset the f vector
        f_c_c = zeros(N,1);
        F_c = sparse(N,N);
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
            f_c_e = @(y) - c(y);
            % enumerate the basisfunctions on interval k.
            enum([1 2]) = [k k+1];
            for i=1:2
                % evaluate intergrals related to f.
                fmul1 = @(y) basis{i}(y) .* f_c_e(y);
                f_c_c(enum(i) ) = f_c_c(enum(i)) + integral(fmul1, x1,x2);
                for j=1:2
                    fmul3 = @(y) c(y).^2 .* basis{i}(y) .* basis{j}(y); %not sure if right
                    F_c(enum(i),enum(j)) = F_c(enum(i),enum(j)) + integral(fmul3,x1,x2);
                
    
            end
            end
        end
    
    end

elseif splitting == 2
        while counter <10000
        counter = counter + 1;

        subplot(2,2,1)
        plot(x, v_n)
        ylabel("composition u");
        xlabel("x");
        ylim([-1.5,1.5])
        pause(0.00001)

        G = [M, dt.* A; -epsilon.^2.* A, M];
        b = [M*v_n;  f];
        sol = G \ b; 
        v_n = sol(1:N);
        t(counter) = t(counter-1) + dt;

        subplot(2,2,2)
        mass(counter) = sum(v_n)./N;
        plot(t, mass)
        u_n = sol(N+1:end);
        xlabel("t");
        ylabel("mass");

        subplot(2,2,3)
        E(counter) = free_engergy(v_n);
        plot(t, E)
        ylabel("total free energy");
        xlabel("t");
        
    
        %reset the f vector
        f_c_c = zeros(N,1);

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
            f_all = @(y) c(y).^3 - c(y); % for other splitting
            % enumerate the basisfunctions on interval k.
            enum([1 2]) = [k k+1];
            for i=1:2
                % evaluate intergrals related to f.
                fmul11 = @(y) basis{i}(y) .* f_all(y);
                f(enum(i)) = f_c_c(enum(i)) + integral(fmul11, x1,x2);
            end
        end
    
        end

        %fully explicit scheme
elseif splitting == 3
        while counter <10000
            counter = counter + 1;
    
            subplot(2,2,1)
            plot(x, v_n)
            ylabel("composition u");
            xlabel("x");
            ylim([-1.5,1.5])
            pause(0.00001)
    
            G = [M, dt.* A; sparse(N,N), M];
            b = [M*v_n ;  f + epsilon.^2.* A * v_n];
            sol = G \ b; 
            v_n = sol(1:N);
            t(counter) = t(counter-1) + dt;
    
            subplot(2,2,2)
            mass(counter) = sum(v_n)./N;
            plot(t, mass)
            u_n = sol(N+1:end);
            xlabel("t");
            ylabel("mass");
    
            subplot(2,2,3)
            E(counter) = free_engergy(v_n);
            plot(t, E)
            ylabel("total free energy");
            xlabel("t");
            
        
            %reset the f vector
            f_c_c = zeros(N,1);
    
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
                f_all = @(y) c(y).^3 - c(y); % for other splitting
                % enumerate the basisfunctions on interval k.
                enum([1 2]) = [k k+1];
                for i=1:2
                    % evaluate intergrals related to f.
                    fmul11 = @(y) basis{i}(y) .* f_all(y);
                    f(enum(i)) = f_c_c(enum(i)) + integral(fmul11, x1,x2);
                end
            end
    
        end

end


     




% b = vhat(2:(N-1),1);
% u(1,1) = 0;
% u(N,1) = 0;
% u(2:(N-1),1) = A\b