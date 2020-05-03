clf; clc; close all; clear all;
pkg load statistics;

%
% Step 1: Process Data
%

% Step 1.1: Load Population Data

population_spain   = 46450000; % Statista.com (2018)
N                  = population_spain;

% Step 1.2: Load Infected Cases

[infected1,infected2] = textread('Spain-Infected.csv','%s;%f','delimiter',';');
infected2             = infected2(end:-1:1);
dates                 = (1:1:length(infected2))';

% Step 1.3: Load Recovered Cases

[recovered1,recovered2] = textread('Spain-Recovered.csv','%s;%f','delimiter',';');
recovered2              = recovered2(end:-1:1);

% Step 1.4: Load Dead Cases

[dead1, dead2] = textread('Spain-Dead.csv','%s;%f','delimiter',';');
dead2          = dead2(end:-1:1);

% Step 1.5: Process Data according to (7) from Preprint

R              = recovered2 + dead2;
I              = infected2 - R;
S              = N - I - R;

%
% Step 2: Solve Linear System
%

% Step 2.1: Build Solution Vectors Alpha and Beta

Alpha = zeros(length(S),1);
Beta  = zeros(length(R),1);
Alpha(1) = 0.0;
Beta(1)  = 0.0;

% Step 2.2: Solve 

for j = 1:1:(length(S)-1)
  if S(j + 1) != 0 & I(j + 1) != 0
    Alpha(j + 1) = N*(S(j) - S(j + 1))/(S(j + 1)*I(j + 1));
    Beta(j + 1)  = (R(j + 1) - R(j))/(I(j + 1));
  elseif S(j + 1) == 0 & I(j + 1) != 0
    Alpha(j + 1) = 0;
    Beta(j + 1)  = (R(j + 1) - R(j))/(I(j + 1));
  else
    Alpha(j + 1) = 0;
    Beta(j + 1)  = 0;
  endif
endfor

%
% Step 3: Statistical Analysis of Beta + Expontential Model for Alpha
%

% Step 3.1: Statistical Analysis of Beta
Stat_Start  = 35;
Beta_Mean   = sum(Beta(Stat_Start:1:end))/(length(Beta(Stat_Start:1:end))-1);
Beta_Median = median(Beta(Stat_Start:1:end));
Beta_SD     = ((1/(length(Beta(Stat_Start:1:end))-1))*sum((Beta(Stat_Start:1:end)-Beta_Mean)).^2)^0.5;

Beta_Mean_Vec   = Beta_Mean*ones(length(dates),1);
Beta_Median_Vec = Beta_Median*ones(length(dates),1);

% Step 3.2: Exponential Model for Alpha

dates_mod   = dates(35:1:end);
Alpha_mod   = Alpha(35:1:end);

x_Alpha_SSD = fminsearch (@(x) ( sum( abs((x(1)*exp(-x(2)*dates_mod) - Alpha_mod) )) ),[0.8;0.2]);
Alpha_SSD   = x_Alpha_SSD(1)*exp(-x_Alpha_SSD(2)*dates);

% Step 3.3: Weibull Model for Beta

%x_Beta_SSD  = fminsearch (@(x) ( sum( (abs( (x(1)/x(2))*((dates/x(2)).^(x(1)-1)).*exp(-((dates/x(2)).^(x(1)-1))) - Beta) ).^2) ),[1.0;1.0]);
%Beta_SSD    = (x_Beta_SSD(1)/x_Beta_SSD(2))*((dates/x_Beta_SSD(2)).^(x_Beta_SSD(1)-1)).*exp(-(-(dates/x_Beta_SSD(2)).^(x_Beta_SSD(1)-1)));

% Step 3.4: Gaussian Model for Beta
x_Beta_SSD   = fminsearch (@(x) ( sum( abs( x(1)*exp(-(((dates-x(2))/x(3)).^2)) - Beta ) ) ),[0.25;65.0;20.0]);
Beta_SSD     = x_Beta_SSD(1)*exp(-((dates-x_Beta_SSD(2))/x_Beta_SSD(3)).^2);

%
% Step 4: Plots
%

% Step 4.1: Plot Infected and Recovered + Solutions

figure(1)
hold off
plot(dates, I, 'color', 'black', '-@');
lt01 = title("Spain: People vs. Dates", "fontsize", 14);
xt01 = xlabel("t", "fontsize", 14);
yt01 = ylabel("People", "fontsize", 14);
hold on
plot(dates, R, 'color', 'blue', '-@');
legend({"Infected", "Recovered"}, 'location', 'northwest');

figure(2)
hold off
plot(dates, Alpha, 'color', 'red', '-@');
lt02 = title("Spain: Transfer Rates vs. Dates", "fontsize", 14);
xt02 = xlabel("t", "fontsize", 14);
yt02 = ylabel("Transfer Rates", "fontsize", 14);
hold on
plot(dates, Beta, 'color', 'green', '-@');
ylim([0.0 0.8]);
legend({"Alpha", "Beta"}, 'location', 'northwest');

figure(3)
hold off
plot(dates(35:1:end), Alpha(35:1:end), 'color', 'red', '-@');
lt03 = title("Spain: Alpha vs. Dates", "fontsize", 14);
xt03 = xlabel("t", "fontsize", 14);
yt03 = ylabel("Transfer Rate Alpha", "fontsize", 14);
hold on
plot(dates(35:1:end), Alpha_SSD(35:1:end), 'color', 'magenta', '-@');
ylim([0.0 0.8]);
legend({"Alpha (Data)", "Alpha (Exponential)"}, 'location', 'northwest');

figure(4)
hold off
plot(dates, Beta, 'color', 'green', '-@');
hold on
plot(dates, Beta_Mean_Vec, 'color', 'yellow', '-@');
hold on
plot(dates, Beta_Median_Vec, 'color', 'cyan', '-@');
hold on
plot(dates, Beta_SSD, 'color', 'magenta', '-@');
ylim([0.0 0.4]);
lt04 = title("Spain: Beta vs. Dates", "fontsize", 14);
xt04 = xlabel("t", "fontsize", 14);
yt04 = ylabel("Transfer Rate Beta", "fontsize", 14);
ht04 = legend ({"Beta", "Mean", "Median", "Beta (Gaussian)"}, 'location', 'northwest');

figure(5)
hold off

subplot(2,2,1);
hold off
plot(dates, I, 'color', 'black', '-@');
lt01 = title("Spain: People vs. Dates", "fontsize", 14);
xt01 = xlabel("t", "fontsize", 14);
yt01 = ylabel("People", "fontsize", 14);
hold on
plot(dates, R, 'color', 'blue', '-@');
legend({"Infected", "Recovered"}, 'location', 'northwest');

subplot(2,2,2);
hold off
plot(dates, Alpha, 'color', 'red', '-@');
lt02 = title("Spain: Transfer Rates vs. Dates", "fontsize", 14);
xt02 = xlabel("t", "fontsize", 14);
yt02 = ylabel("Transfer Rates", "fontsize", 14);
hold on
plot(dates, Beta, 'color', 'green', '-@');
ylim([0.0 0.8]);
legend({"Alpha", "Beta"}, 'location', 'northwest');

subplot(2,2,3);
hold off
plot(dates(30:1:end), Alpha(30:1:end), 'color', 'red', '-@');
lt03 = title("Spain: Alpha vs. Dates", "fontsize", 14);
xt03 = xlabel("t", "fontsize", 14);
yt03 = ylabel("Transfer Rate Alpha", "fontsize", 14);
hold on
plot(dates(35:1:end), Alpha_SSD(35:1:end), 'color', 'magenta', '-@');
ylim([0.0 0.8]);
legend({"Alpha (Data)", "Alpha (Exponential)"}, 'location', 'northwest');

subplot(2,2,4);
hold off
plot(dates, Beta, 'color', 'green', '-@');
hold on
plot(dates, Beta_Mean_Vec, 'color', 'yellow', '-@');
hold on
plot(dates, Beta_Median_Vec, 'color', 'cyan', '-@');
hold on
plot(dates, Beta_SSD, 'color', 'magenta', '-@');
ylim([0.0 0.4]);
lt04 = title("Spain: Beta vs. Dates", "fontsize", 14);
xt04 = xlabel("t", "fontsize", 14);
yt04 = ylabel("Transfer Rate Beta", "fontsize", 14);
ht04 = legend ({"Beta", "Mean", "Median", "Beta (Gaussian)"}, 'location', 'northwest');