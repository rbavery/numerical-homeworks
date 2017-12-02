% 1. warm up

%1
% fxa roots: (1.506, 0)
fxa = @(x)exp(x) - x - 3;
fxaroot = NewtonMethod1(fxa, 1,7);
fxaroot2 = NewtonMethod1(fxa, -3,11)
% fxb roots: (1.2208, 0)
fxb = @(x)exp(2.*x) - x.^2 - 10;
fxbroot = NewtonMethod1(fxb, 1,8);
% fxc roots: (1.5567, 0)
fxc = @(x)log(x) + x - 2;
fxcroot = NewtonMethod1(fxc, 1,9);
% fxd roots: (3.2246, 0)
fxd = @(x)x.*log(x) + x - 7;
fxdroot = NewtonMethod1(fxd, 1,10);

%2
%a
z = [100 250 650 1000 1800 2350 3500 5000];
ocf = [5.8166 2.4736 1.4481 1.0013 0.2982 0.4126 0.2739 0.2897];
figure(1)
plot(ocf, z, 'b')
hold on
title('Flux as function of depth')
xlabel('Organic Carbon Flux')
ylabel('Depth')



%b

znought = 100;
b = .01;
[c, cprime, cprimeprime] = costfun(b, z, ocf, .001);
counter = 1;

while abs(cprime) >= 1e-2
    bi(counter) = b
    ci(counter) = c;
    cprimei(counter) = cprime;
    b = b - cprime./cprimeprime;
    
    [c, cprime, cprimeprime] = costfun(b, z, ocf, .001);
    
    
    counter = counter + 1;
     
    end

b
Fz = 5.8166.*(z./100).^-b;
plot(Fz,z, 'r');

cprimeprime > 0


% c

figure(2)
subplot(2, 1, 1);
plot(bi, ci);
title('b guess .01');
hold on
xlabel('b value');
ylabel('cost');

subplot(2, 1, 2);
plot(bi, cprimei);
xlabel('b value');
ylabel('cost derivative');
hold off

%d

b2 = 1.9;
[c2, cprime2, cprimeprime2] = costfun(b2, z, ocf, .001);
counter2 = 1;

while abs(cprime2) >= 1e-2
    bi2(counter2) = b2
    ci2(counter2) = c2;
    cprimei2(counter2) = cprime2;
    b2 = b2 - cprime2./cprimeprime2;
    
    [c2, cprime2, cprimeprime2] = costfun(b2, z, ocf, .001);
    
    
    counter2 = counter2 + 1;
     
    end

figure(3)    
subplot(2, 1, 1);
plot(bi2, ci2);
hold on
xlabel('b value');
ylabel('cost');
title('b guess 1.9');

subplot(2, 1, 2);
plot(bi2, cprimei2);
xlabel('b value');
ylabel('cost derivative');

b2
Fz2 = 5.8166.*(z./100).^-b2;
plot(Fz2,z, 'r');


function [c, cprime, cprimeprime] = costfun(barg, z, ocfobs, db)
    
    % the model
    Fz = @(barg)5.8166.*(z./100).^-barg;
    
    % the cost function
    cost = @(barg)sum((Fz(barg)-ocfobs).^2);
    
    c = cost(barg);
    
    cup = cost(barg+db);
    cdown = cost(barg-db);
    
    cprime = (cup - cdown)./(db*2);
    
    cprimeup = (cup - c)./(db); 
    cprimedown = (c - cdown)./(db);
    
    cprimeprime = (cprimeup - cprimedown)./(db);
    
end    

function x = NewtonMethod1(fx, xguess, i)
    y = fx(xguess);
    counter = 1
    dx = .1;
    while  abs(y) >= 1e-2
        xi(counter) = xguess;
        yi(counter) = y;
        % centered difference
        dfdx = (fx(xguess+dx)-fx(xguess-dx))./(2*dx);
        xguess = xguess - fx(xguess)./dfdx;
        y = fx(xguess);
        counter = counter + 1;
    end
    figure(i);
    plot(xi, yi)
    xlabel('x')
    ylabel('y')
    x = xguess;
end