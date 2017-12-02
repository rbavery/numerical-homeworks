% open data
M = importdata('API_SP.POP.TOTL_DS2_en_csv_v3.csv');
data = M.data;
headers = M.textdata;

% the year is the first row
% the rest of the rows are the population data
year = data(1,:);

% find the relevant countries (the index of the row for each country)
ilow = find(strcmp(M.textdata,'Low income'));
ilowmid = find(strcmp(M.textdata,'Lower middle income'));
ihigh = find(strcmp(M.textdata,'High income'));
ijapan = find(strcmp(M.textdata,'Japan'));
ichina = find(strcmp(M.textdata,'China'));
ius = find(strcmp(M.textdata,'United States'));

% now you can make the plots, answer the questions etc.
% 1a
figure(1)
x = [1960:1:2016];
dlow = data(ilow, :)./1000000;
dlowmid = data(ilowmid, :)./1000000;
dhigh = data(ihigh, :)./1000000;
semilogy(x, dlow, 'b');
hold on
% e
linearparams = polyfit(x, log(dlow), 1);
linearfit_low = exp(linearparams(1).*x+linearparams(2));
semilogy(x, linearfit_low, 'm');

semilogy(x, dlowmid, 'r');
semilogy(x, dhigh, 'g');
title('Population Levels by Year');
legend('Low Income','Linear Fit to Low Pop', 'Lower Middle Income','High Income');
xlabel('year');
ylabel('population in millions (log scale)');
hold off

% 1b growth rates
figure(2)

dt = 1;
change_dlow = dlow(2:end) - dlow(1:end-1);
change_dlowmid = dlowmid(2:end) - dlowmid(1:end-1);
change_dhigh = dhigh(2:end) - dhigh(1:end-1);
semilogy(x(2:end), change_dlow, 'b');
hold on
semilogy(x(2:end), change_dlowmid, 'r');
semilogy(x(2:end), change_dhigh, 'g');
title('Population Growth by Year');
legend('Low Income','Lower Middle Income','High Income');
xlabel('year');
ylabel('population in millions (log scale)');
hold off


% 1c
% low income is growing exponentially
% can tell because trend is linear in the log domain

% 1d
% 1985ish

%1e
% linearparams = polyfit(x, log(dlow), 1);
% linearfit_low = exp(linearparams(1).*x+linearparams(2));
% a = .0261 by logging both sides of exponential growth model
% then fitting model in matlab and getting slope and intercept
% got P nought from the data
Pnought = dlow(1)

%f
Pt = 154670163.*exp(.0261.*x);

%g
Pt2050 = 154670163.*exp(.0261.*(2050-1960))./1000000

%2a

figure(3)
djapan = data(ijapan, :)./1000000;
plot(x, djapan, 'b');
hold on

% 2b
p = polyfit(x, djapan, 2);
%equation below
p(1)
p(2)
p(3)
seconddegfit = (p(1).*x.^2)+p(2).*x+p(3);
plot(x, seconddegfit, 'm');
title('Population Levels of Japan by Year');
legend('Japan Population','Linear Fit', 'Location', 'northwest');
xlabel('year');
ylabel('population in millions');

%b and c
seconddegfit = p(1)*2050^2+p(2)*2050+p(3)
dx = 2050 - 1960;
taylorseries2 = p(1)*1960^2+p(2)*1960+p(3)+(p(1)*2*1960+p(2))*dx
taylorseries3 = p(1)*1960^2+p(2)*1960+p(3)+(p(1)*2*1960+p(2))*dx+(p(1)*dx^2)

dx = 2016 - 1960;
taylorseries2016 = p(1)*1960^2+p(2)*1960+p(3)+(p(1)*2*1960+p(2))*dx+(p(1)*dx^2)


























