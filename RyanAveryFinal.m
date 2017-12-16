% 1
x = [0:1:400];
expr = ((x + 1).*(2*x+3))./((x + 2).*(3*x+4));
limitvalue(1:1:401) = .666666666666666;

figure(1)
plot(x, expr);
hold on
plot(x, limitvalue);
xlabel('x')
ylabel('expression as x approaches infinity')
title('Limit of Problem 1 Expression')
legend('expression', 'the limit as x approaches infinity');

%2
dt = .1;
t = [0:dt:10];
x = exp(t./2).*sin(2.*t);
y = exp(t./5).*cos(2.*t);

figure(2)
plot(x , y);
xlabel('x')
ylabel('y')
title('2a Position in Time Domain')
legend('position');

unum = (x(2:end)-x(1:end-1))./dt;
vnum = (y(2:end)-y(1:end-1))./dt;
u = exp(t./2).*(.5*sin(2*t)+2*cos(2.*t));
v = exp(t./5).*((1/5)*cos(2*t)-2*sin(2.*t));
figure(3)
subplot(2, 1, 1);
plot(unum, vnum);
title('2d Derivatives with Numerical Method (top) and Analytical Method (bottom)');
xlabel('u');
ylabel('v');
hold on

subplot(2, 1, 2);
plot(u, v);
xlabel('u');
ylabel('v');
hold off

%3a
dt = 100;
t = [0:dt:1000];
tmids = (t(2:end)+t(1:end-1))/2;
Ag = [535 390 235 166 108 65 47 28 22 12 8];

params = polyfit(t, log(Ag), 1);
lnA0 = params(2);
k = params(1);
Ag = exp(lnA0+k*t);

%3b half life
thalf = (log(Ag(1)/2)-lnA0)/k;
Ag2000 = exp(lnA0+k*2000);

% 3cd
dAgdta = exp(lnA0)*k*exp(k*tmids);
dAgdt = diff(Ag)/dt;

figure(4)
plot(t,Ag)
xlabel('Year')
ylabel('Grams')
legend('Exponential Model')
title('3a Radioactive Decay of Element') 

figure(5)
plot(tmids,dAgdta, '-.b')
hold on
plot(tmids,dAgdta,'r*');
xlabel('Year')
ylabel('Change in Grams (by 100 years)')
legend('analytic', 'numerical')
title('3d Derivative of the Radioactive Decay of Element') 

%3e)
% analytically for 0 to 1000
intAga = -134857*exp(-.0042*1000)+134857;
% numerically
intAgn = sum(dt*Ag);

%4c)
dx = .01;
x = [-5:dx:5];
fx = 4*x + 4 - exp(x/2)*4;
Fx = 1.405+.7026.*( x-1 )-.8244.*( x-1 ).^2;

figure(6)
plot(x,fx, 'b')
hold on
plot(x,Fx,'r');
xlabel('x')
ylabel('y')
legend('f(x)', 'F(x) (approximation)')
title('4c Comparison of f(x) and Taylor Series approx')

% 4d)

x = 1.5;
f = @(x)(4*x + 4 - exp(x/2)*4);

y = f(x);
itc = 0;
while abs(y) > 1e-3
    
    y = f(x);
    fup = f(x+dx);
    fdown = f(x-dx);
    fprime = (fup-fdown)./(2*dx);
    
    x = x - y./fprime;
    itc = itc + 1;
    xi(itc) = x;
    yi(itc) = y;
end
itc
xi
yi
x
    
%5a
dx = .1;
x = [-6:dx:6];
y = [-6:dx:6];
[X, Y] = meshgrid(x, y);
z = 2.*sin(X./2)+cos(Y)+exp(X./5);
%5b
dzdx = @(x)cos(x/2)+.2*exp(x/5);
dzdy = @(y)-sin(y);
x = -4;
y = 2;
zx=dzdx(x);
zy=dzdy(y);

while abs(zx) > 1e-3 || abs(zy) > 1e-3
    zx=dzdx(x);
    zxup = dzdx(x+dx);
    zxdown = dzdx(x-dx);
    zxprime = (zxup-zxdown)./(2*dx);
    x = x - zx./zxprime;
    
    zy=dzdy(y);
    zyup = dzdy(y+dx);
    zydown = dzdy(y-dx);
    zyprime = (zyup-zydown)./(2*dx);
    y = y - zy./zyprime;
end

zx;
zy;
x;
y;
figure(7)
surf(X, Y, z);
xlabel('x')
ylabel('y')
zlabel('z')
title('Surf plot')
colorbar()

% 5d
dzdX = dzdx(X);
dzdY = dzdy(Y);
figure(8)
subplot(2, 1, 1)
contourf(X, Y, dzdX);
xlabel('x')
ylabel('y')
title('5d Partial Derivative By X')
subplot(2 , 1, 2)
contourf(X, Y, dzdY);
xlabel('x')
ylabel('y')
title('5d Partial Derivative By Y')
colorbar()

%6c
dt = 0.1;
t = [0:dt:10];

k = .5;
ym = 10;
y0 = 1;
dydt = @(y)k*(ym-y);
%Euler Forward
ef = 0*t; % placeholder
ef(1) = y0;

for i = 2:length(t)
    ef(i) = ef(i-1) + dt*dydt(ef(i-1));
end

%Euler Backward
bf = 0*t;
bf(1) = y0;

for i = 2:length(t)
    % rearrange terms for backward
    % bf(i) = bf(i-1) + dt*.5*(10-bf(i))
    % bf(i) + .5*bf(i) = bf(i-1) + dt*5
    bf(i) = (bf(i-1) + dt*5)./1.5;
end

% Midpoint
mf = 0*t;
mf(1) = y0;
for i = 2:length(t)
    % Euler forward 1/2 step
    mf(i) = mf(i-1) + dt./2*dydt(mf(i-1));
    % full step with midpoint
    mf(i) = mf(i-1) + dt*dydt(mf(i));
end

f = @(t)-9.*exp(-k*t)+ym; % analytic solution

figure(9)
plot(t, ef, '-*')
hold on
plot(t, bf, '-x')
hold on
plot(t, mf, '-o')
hold on
plot(t, f(t), '-.')
legend('forward', 'backward', 'midpoint', 'solved')
xlabel('x')
ylabel('y')
title('6c Solutions to Differential Equations')

%7a
% inputs
a = .2;
p = .1;
q = .1;
b = .3;
r = .2;
s = .1;
c = 5;
d = 10;

% linear terms
A = [.2  0; 0 .3];

% constants
r = [5;10];

%initial guess
x0 = 10;
y0 = 10;

% Using Euler Forward
v0 = [x0;y0]; % initial condition

% time domain
dt = 0.01;
t = [0:dt:6];

v1 = zeros(2, length(t)); % placeholder for solution
v1(:,1) = v0;


% Euler Forward
for i = 2:length(t)
    x = v1(1, i-1);
    y = v1(2, i-1);
    q = [-.1.*x.^2-.1.*x.*y;-.2.*x.*y-.1.*y.^2]; % non-linear terms
    v1(:,i) = v1(:,i-1) + dt*(A*v1(:,i-1) + q + r); 
end
x
y
figure(10)
plot(t, v1)
xlabel('time')
ylabel('population')
title('7b Steady state of populations x and y')
legend('x', 'y')

% intitial conditions
v = [10; 10]; % v = [x0; y0]
% first F value with guess
F = Jacobianfunc(v, A, r); 

while norm(F) > 1e-3
    
    %evaluate function and Jacobian
    [F, J] = Jacobianfunc(v, A, r);
    
    % update v with Newton step, approx
    v = v - (J\F);

    % evaluate function and Jacobian
    F = Jacobianfunc(v, A, r);
    
end

v

function [F, J] = Jacobianfunc(v, A, r)

x = v(1);
y = v(2);

q = [-.1.*x.^2-.1.*x.*y;-.2.*x.*y-.1.*y.^2]; % non-linear terms

% Evaluate the function

F = A*v + q + r;

% Evaluate Jacobian J = A + dq/dv
dqdv = [-.2.*x-.1.*y, -.1.*x; -.2.*y, -.2.*x-.2.*y]; 
J = A + dqdv;
end
    
    
    
    
    
    
    















