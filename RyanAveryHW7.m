% Question 1
% dN/dt = qp*P - kr*N/(N+kN) + S
% dP/dt = kr*N/(N+kN) - a*P*C - qp*P
% dC/dt = a*P*C - qc*C - w*C

% linear terms
% Q R S 
A = [1  0; 0 1];

% constants
r = [S;0;0];

%initial guess
N0 = 1;
P0 = 1;
C0 = 1;

v = [N0; P0; C0]; % v = [x0; y0]

% iterate to solution using Newton's method
F = myfun9(v, A, r); 