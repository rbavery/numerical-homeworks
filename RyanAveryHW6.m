dt = 1/365; % 1 day time step
t = [0:dt:1000]; % 1 years
area = 1000 *1000;
%precip
Pmax = 1; % meter per year
P = (Pmax/2) * (cos(2*pi*t)+1); % meters cubed of precip over 1km area

% figure(1)
% plot(t, P, '-b')
% title('Precip')

%temp

Tmax = 25; %deg celsius
Tmin = 10;
T = ((Tmax-Tmin)/2)*(sin(2*pi*t-pi/2)+1)+Tmin;

% %plot
% figure(2)
% plot(t,T)
% title('Temp')

%runoff = q_s
ks = 10; % no units
Ph = 1; % meter per year
qs = (P./(1+ks*exp(-P./Ph))); %runoff (m/yr) make sure to multiple by area to get cubic meters per year

% figure(3)
% plot(P, qs./P)
% title('Fraction of Precip that goes to runoff')

%1 see doc

%2
kt = .001; % meters per year*celsius^2
Pmax = 1; % meters cubed per year
Tmax = 25; %deg celsius
Tmin = 10;

E = (kt*.5*T.^2); % meters cubed of evap over 1km area

% figure(2)
% plot(t(1:365), P(1:365), 'r')
% hold on
% plot(t(1:365), E(1:365), 'b')
% title('Precipitation and Evaporation over 1 Year over 1km aquifer')
% xlabel('Time (in Proportion of Year by Days)')
% ylabel('Meters of precip')

totalrainfall = sum(P*dt); % .5027 meters per one year

%3 
qn = 0; % no groundwater pumping
kb = .01; % per year

ks = 10; % no units
Ph = 1; % meter cubed per year


% initial condition
S0 = 10^6; % meters cubed
% solve with Euler Forward

Sf = 0*t; % placeholder

Sf(1) = S0;

f = @(Sf, P, E, qs)P*area-E*area-qs*area-(kb.*Sf);

for i = 2:length(t)
    Sf(i) = Sf(i-1) + dt.*f(Sf(i-1), P(i-1), E(i-1), qs(i-1));
end

steadyyear = 800;
startsteadyi = 800*365;
endfiftyyearsi = startsteadyi+50*365;

rangestorage = range(Sf(startsteadyi:endfiftyyearsi))
meandepth = (mean((Sf(startsteadyi:endfiftyyearsi)))/1000)/1000

figure(3)
plot(t(startsteadyi:endfiftyyearsi), Sf(startsteadyi:endfiftyyearsi))
title('Water Storage at Steady State')
xlabel('Time (in Proportion of Year by Days)')
ylabel('Meters Cubed')

% 4
qn = 100000; % m cubed per year

% initial condition
S0 = Sf(startsteadyi); % meters cubed
% solve with Euler Forward

Sf = 0*t; % placeholder

Sf(1) = S0;

f = @(Sf, P, E, qs, qn)P*area-E*area-qs*area-(kb.*Sf)-qn;

for i = 2:length(t(startsteadyi:end))
    Sf(i) = Sf(i-1) + dt.*f(Sf(i-1), P(i-1), E(i-1), qs(i-1), qn);
end

figure(4)
plot(t(1:365*50), Sf(1:365*50))
hold on
plot(t(365*50:end), Sf(365*50:end))
title('Water Storage Steady State with Pumping for Fifty Years')
xlabel('Time (in Proportion of Year by Days)')
ylabel('Meters Cubed')


% 5

qn0 = 0; % initial pumping m cubed per year
qn = 0*t;
qn(1) = qn0;

for i = 2:length(t)
    qn(i) = qn(i-1) + dt.*10000;
end

Sf = 0*t; % placeholder

Sf(1) = S0;

for i = 2:length(t)
    Sf(i) = Sf(i-1) + dt.*f(Sf(i-1), P(i-1), E(i-1), qs(i-1), qn(i-1));
end

figure(5)
plot(t(startsteadyi:endfiftyyearsi), Sf(startsteadyi:endfiftyyearsi))
title('Water Storage when Pumping Rate Grows with Pop')
xlabel('Time (in Proportion of Year by Days)')
ylabel('Meters Cubed')










