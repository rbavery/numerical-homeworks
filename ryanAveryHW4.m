% open U.S. survival curve data
fid = fopen('us_survival_data_2013.txt');
f = textscan(fid,'%s%s%f%f%f%f%f%f%f%f','HeaderLines',2,...
    'Delimiter',',');
fclose(fid);

% numbers for 2013
la_f = f{3}; % survivorship curve for females in 2013 (number of people
             % surviving to the beginning of that age interval)
la_m = f{4}; % survivorship curve for males in 2013


% 1a Calc life expect as a func of age for males and females in 2013. Plot e as function of a for both.

a = [0 1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105];
da = a(2:end)-a(1:end-1);
[female_a, female_e] = lifeexpect(da, la_f, a);
[male_a, male_e]  = lifeexpect(da, la_m, a);

figure(1);
plot(female_a, female_e, '.-r');
title('Life Expectancy between Females (red) and Males (green)');
hold on
plot(male_a, male_e, '.-g');
xlabel('Age');
ylabel('Life Expectancy');
legend('Female','Male');

% 1b

femalediff = abs(female_e - transpose(female_a));
findex = femalediff < 5;
% 37.5 for females
fresult = female_a(findex)

% 37.5 for males
malediff = abs(male_e - transpose(male_a));
mindex = malediff < 1;
mresult = male_a(findex)

% 2a

D_prob_f = probdeath(da, la_f, a);
D_prob_m = probdeath(da, la_m, a);

figure(2)
semilogy(female_a, D_prob_f, '.-r');
hold on
semilogy(female_a, D_prob_m, '.-g');
xlabel('Age');
ylabel('Probability of Death');
title('Probability of Death for Females (red) and Males (green)');
legend('Female', 'Male');
at = transpose(a(1:end-1))
size(at)
size(D_prob_f)
size(D_prob_m)
fem_death_area = sum(D_prob_f*5) % 1.0214 prob of death occuring between start (2.5) and end point (102.5)
male_death_area = sum(D_prob_m*5) % 1.0259. Can't use step array da for some reason, large values >1

% highest prob of death for males at 97.5 and females at 100
% found by printing probs to console and matching indexes of max prob with xmids array
D_prob_f;
D_prob_m;
female_a(end-3);



%3a
la_fmids = .5 * (la_f(1:end-1)+la_f(2:end));
la_mmids = .5 * (la_m(1:end-1)+la_m(2:end));
u_f = D_prob_f./la_fmids;
u_m = D_prob_m./la_mmids;

figure(3)
subplot(2,1,1)
semilogy(female_a, u_f, 'r');
hold on
semilogy(male_a, u_m, 'g');
xlabel('Age');
ylabel('Force of Mortality');
title('Force of Mortality for Females (red) and Males (green)');
legend('Female', 'Male')
subplot(2,1,2);
reldiff = (u_m - u_f)./u_f;
semilogy(female_a, reldiff, '--');
legend('Relative Difference')

%3b
maxindex = find(reldiff == max(reldiff));
female_a(maxindex) % at maximum 102.5?
% 
% 
%4a
da = transpose(da);
C_f = exp(cumsum(da.*(-1)./(female_e)));
C_m = exp(cumsum(da.*(-1)./(male_e)));

figure(4);
plot(female_a, C_f, '.-r');
title('Proportion of Population Age a and Older, Females and Males');
hold on
plot(male_a, C_m, '.-g');
xlabel('Age');
ylabel('Proportion of Pop');
legend('Female','Male');

%b
%found median in janky way with print to console
C_f
medC_f = C_f(9);
female_a(9); % 

C_m
medC_m = C_m(8);
male_a(8); %

%c index into C_* with same index of age closest to 59.5
retindex = 14
C_f(retindex)
C_m(retindex)

function answer = probdeath(dx, yarr, xarr)
    
    dx = transpose(dx);
    D = -1./yarr(1).*(yarr(2:end)-yarr(1:end-1))./dx;
    
    
    answer = D;
end

function [xbmids, answer] = lifeexpect(dx, yarr, xarr)
    
    nxb = length(xarr);
    ny = length(yarr);
    
    xmidpoints = @(x, nx).5 * (x(1:nx-1)+x(2:nx));
    xbmids = xmidpoints(xarr, nxb);
    
    ymidpoints = @(y, n).5 * (y(1:n-1)+y(2:n));
    ybmids = ymidpoints(yarr, ny);
    
    dx = transpose(dx);
    
    y_int_0_105 = sum(ybmids.*dx)./ybmids(1);
    y_int_0_a = cumsum(ybmids.*dx)./ybmids(1);
    y_int_a_105 = y_int_0_105 - y_int_0_a;
     
    answer = y_int_a_105;
    
end