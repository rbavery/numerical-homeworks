% Week 6 Section

% euler methods

dt = 0.1;
t = [0:dt:10];

% differential equation to solve
f = @(t,y).5.*exp(-t./2)-.5.*y;
solvedf3 = @(t).5.*t.*exp(-t./2);

asolvedf3 = solvedf3(t);
% initial condition
y0 = .5;

% solve with Euler Forward

yf = 0*t; % placeholder


yf(1) = y0;


for i = 2:length(t)
    yf(i) = yf(i-1) + dt*f(t(i-1),yf(i-1));
end

% solve with Euler Backward
yb = 0*t; % placeholder

yb(1) = y0; % initial condition


for i = 2:length(t)
    yb(i) = yb(i-1)./(1-dt.*cos(t(i)));
end

% plot
figure(1)
plot(t, yf, '-ob')
hold on
plot(t, yb, '-xr')
plot(t, asolvedf3, '-g')
legend('forward', 'backward', 'solved')
