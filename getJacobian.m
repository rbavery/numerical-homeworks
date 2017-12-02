function [F, J] = getJacobian(v, A, r)

x = v(1);
y = v(2);

q = [-x.^2-x*y; -x*y-y.^2]; % non-linear terms

% Evaluate the function

F = A*v + q + r;

% non-linear terms
%- Px^2 - Qxy
%- Rxy, - Sy^2

% Evaluate Jacobian J = A + dq/dv
dqdv = [-2*x-y, -y; -y, -x-2*y]; 
J = A + dqdv;
end