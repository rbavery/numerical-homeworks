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

figure(2)
plot(t(1:365), P(1:365), 'r')
hold on
plot(t(1:365), E(1:365), 'b')
title('Precipitation and Evaporation over 1 Year over 1km aquifer')
xlabel('Time (in Proportion of Year by Days)')
ylabel('Meters of precip')

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

f = @(Sf, P, E, qs, qn)P*area-E*area-qs*area-(kb.*Sf)-qn;

for i = 2:length(t)
    Sf(i) = Sf(i-1) + dt.*f(Sf(i-1), P(i-1), E(i-1), qs(i-1), qn);
end

steadyyear = 800;
startsteadyi = 800*365;
endfiveyearsi = startsteadyi+5*365;
steadystorage = Sf(startsteadyi);

rangestorage = range(Sf(startsteadyi:endfiveyearsi))
meandepth = (mean((Sf(startsteadyi:endfiveyearsi)))/1000)/1000

figure(3)
plot(t(startsteadyi:endfiveyearsi), Sf(startsteadyi:endfiveyearsi))
title('Water Storage at Steady State')
xlabel('Time (in Proportion of Year by Days)')
ylabel('Meters Cubed')

% 4
qn = 100000/365; % m cubed per year

% initial condition

% solve with Euler Forward

Sf4 = 0*t; % placeholder
S0 = steadystorage; % meters cubed
Sf4(1) = S0;

f4 = @(Sf4, P, E, qs, qn)P*area-E*area-qs*area-(kb.*Sf4)-qn;

for i = 2:length(t)
    Sf4(i) = Sf4(i-1) + dt.*f4(Sf4(i-1), P(i-1), E(i-1), qs(i-1), qn);
end

endfiftyyearsi = startsteadyi+50*365;
endcpumping = Sf4(50*365+1);
depthdiff = ((Sf4(startsteadyi)-Sf4(endfiveyearsi))/1000)/1000
figure(4)
plot(t(startsteadyi:endfiftyyearsi), Sf4(1:50*365+1))
hold on
title('Water Storage with Pumping')
xlabel('Time (in Proportion of Year by Days)')
ylabel('Meters Cubed')


% 5

qn0 = 100000; % initial pumping m cubed per year^2
qn = 0*t;
qn(1) = qn0;

qn = qn0 + t*10000;

Sf5 = 0*t; % placeholder

Sf5(1) = endcpumping;

for i = 2:length(t)
    Sf5(i) = Sf5(i-1) + dt.*f4(Sf5(i-1), P(i-1), E(i-1), qs(i-1), qn(i-1));
end

popi = endfiftyyearsi+1 + 50*365;
plot(t(endfiftyyearsi+1:popi), Sf5(1:(1+50*365)))


endpumppop = Sf5(1+50*365);

%6
kn = .0135;
Sf6(1) = endpumppop;

f6 = @(Sf4, P, E, qs)P*area-E*area-qs*area-(kb.*Sf4)-kn*Sf4;

for i = 2:length(t)
    Sf6(i) = Sf6(i-1) + dt.*f6(Sf6(i-1), P(i-1), E(i-1), qs(i-1));
end

regi = popi+ 1 + 50*365;
plot(t(popi+1:regi), Sf6(1:(1+50*365)))
legend('Constant Pumping for Fifty Years', 'Pumping Rate Grows with Pop', 'Pumping Rate Scales with Water Storage')


















