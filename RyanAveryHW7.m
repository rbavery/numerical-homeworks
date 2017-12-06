% Question 1
% dN/dt = qp*P - kr*N/(N+kN) + S
% dP/dt = kr*N/(N+kN) - a*P*C - qp*P
% dC/dt = a*P*C - qc*C - w*C

% 2 inputs
q_p = .5;
q_c = 1/5;
k_R = 2;
k_N = 1/10;
a = 1/3;
S = 1/20;
w = 1/5;

% linear terms
A = [0  q_p q_c; 0 -q_p 0; 0 0 (-q_c-w)];

% constants
r = [S;0;0];

%initial guess
N0 = 1;
P0 = 1;
C0 = 1;

% Using Euler Forward
v0 = [N0;P0;C0]; % initial condition

% time domain
dt = 0.01;
t = [0:dt:200];

v1 = zeros(3, length(t)); % placeholder for solution
v1(:,1) = v0;


% Euler Forward
for i = 2:length(t)
    N = v1(1, i-1);
    P = v1(2, i-1);
    C = v1(3, i-1);
    
    R = k_R.*N./(N+k_N);
    G = a.*P.*C;
    Gi(i) = (G);
    q = [-R;R-G;G]; % non-linear terms
    v1(:,i) = v1(:,i-1) + dt*(A*v1(:,i-1) + q + r); 
end

figure(1)
plot(t, v1)
xlabel('time')
ylabel('units of N')
title('Steady state of N, P, and C')
legend('N', 'P', 'C')

figure(2)
plot(t, Gi)
xlabel('time')
ylabel('units of N')
title('Steady state of G')
legend('G')

%3
% intitial conditions
v = [N0; P0; C0]; % v = [x0; y0]
s = .01;
% first F value with guess
F = getJacobian(v, A, r, k_R, k_N, a); 

%update with Newton Step
itc = 0;
while norm(F) > 1e-3
    
    %evaluate function and Jacobian
    [F, J] = getJacobian(v, A, r, k_R, k_N, a);
    
    % update v with Newton step, approx
    % do not do sJ\Fn in matlab, do s(J\Fn) for proper order
    v = v - s*(J\F);
    
    %update iteration counter
    itc = itc + 1;
    
    % collect values
    Ni(itc) = v(1);
    Pi(itc) = v(2);
    Ci(itc) = v(3);
    normi(itc) = norm(F);
    itci(itc) = itc;
    % evaluate function and Jacobian
    F = getJacobian(v, A, r, k_R, k_N, a);
    
end

figure(3)
plot(itci, normi, '*')
set(gca, 'YGrid', 'on')
set(gca, 'XGrid', 'on')
title('Norm of F, Tolerance = 1e-3')
xlabel('Newton Iteration')
ylabel('norm of Function')

NPC = [Ni; Pi; Ci];
figure(4)
plot(itci', NPC)
set(gca, 'YGrid', 'on')
set(gca, 'XGrid', 'on')
legend('N', 'P', 'C')
title('Np, P, C steady state estimation, Tolerance = 1e-3')
xlabel('Newton Iteration')
ylabel('units of N')
