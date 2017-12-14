% x = [0:1:400];
% exp = ((x + 1).*(2*x+3))./((x + 2).*(3*x+4));
% limitvalue(1:1:401) = .666666666666666;
% figure(1)
% plot(x, exp);
% hold on
% plot(x, limitvalue);
% xlabel('x')
% ylabel('expression as x approaches infinity')
% title('Limit of Problem 1 Expression')
% legend('expression', 'the limit as x approaches infinity');
% dt = .1;
% t = [0:dt:10];
% x = exp(t./2).*sin(2.*t);
% y = exp(t./5).*cos(2.*t);

% figure(2)
% plot(x , y);
% xlabel('x')
% ylabel('y')
% title('Position in Time Domain')
% legend('position');

% unum = (x(2:end)-x(1:end-1))./dt;
% vnum = (y(2:end)-y(1:end-1))./dt;
% u = exp(t./2).*(.5*sin(2*t)+2*cos(2.*t));
% v = exp(t./5).*((1/5)*cos(2*t)-2*sin(2.*t));
% figure(3)
% subplot(2, 1, 1);
% plot(unum, vnum);
% title('Derivatives with Numerical Method (top) and Analytical Method (bottom)');
% xlabel('u');
% ylabel('v');
% hold on
% 
% subplot(2, 1, 2);
% plot(u, v);
% xlabel('u');
% ylabel('v');
% hold off

%3a
t = [0:100:1000];
Ag = [535 390 235 166 108 65 47 28 22 12 8];

params = polyfit(t, log(Ag), 1);
lnA0 = params(2)
k = params(1)
Ag = exp(lnA0+k*t);

%3b half life
thalf = (log(Ag(1)/2)-lnA0)/k
Ag2000 = exp(lnA0+k*2000)

figure(4)
plot(t,Ag)
xlabel('Year')
ylabel('Grams')
legend('Exponential Model')
title('Radioactive Decay of Element') 








