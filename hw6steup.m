dt = 1/365; % 1 day time step
t = [0:dt:10]; % 10 years

%precip
Pmax = 1 % meter per year
P = (Pmax/2) * (cos(2*pi*t)+1);

figure(1)
plot(t, P, '-b')
title('Precip')

%temp

Tmax = 25; %deg celsius
Tmin = 10;
T = ((Tmax-Tmin)/2)*(sin(2*pi*t-pi/2)+1)+Tmin;

%plot
figure(2)
plot(t,T)
title('Temp')

%runoff = q_s
ks = 10; % no units
Ph = 1 % meter per year
qs = P./(1+ks*exp(-P./Ph)); %runoff (m/yr) make sure to multiple by area to get cubic meters per year

figure(3)
plot(P, qs./P)
title('Fraction of Precip that goes to runoff')