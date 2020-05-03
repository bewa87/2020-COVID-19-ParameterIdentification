clf; clc; close all; clear all;
pkg load statistics;

%
% Step 1: Process Data
%

% Step 1.1: Load Population Data

population_hubei = 59170000; % Statista.com (2018)
N                = population_hubei;

% Step 1.2: Load Infected Cases

[infected1,infected2] = textread('China-Hubei-Infected.csv','%s;%f','delimiter',';');
infected2             = infected2(end:-1:1);
dates                 = (1:1:length(infected2))';

% Step 1.3: Load Recovered Cases

[recovered1,recovered2] = textread('China-Hubei-Recovered.csv','%s;%f','delimiter',';');
recovered2              = recovered2(end:-1:1);

% Step 1.4: Load Dead Cases

[dead1, dead2] = textread('China-Hubei-Dead.csv','%s;%f','delimiter',';');
dead2          = dead2(end:-1:1);

% Step 1.5: Process Data according to (7) from Preprint

R              = recovered2 + dead2;
I              = infected2 - R;
S              = N - I - R;

% Step 2: Plot Unprocessed Data

figure(1)
hold off
plot(dates, infected2, 'color', 'black', 'x');
lt01 = title("Hubei (China): People vs. Dates (Unprocessed)", "fontsize", 14);
xt01 = xlabel("t", "fontsize", 14);
yt01 = ylabel("People", "fontsize", 14);
hold on
plot(dates, recovered2, 'color', 'blue', 'o');
hold on
plot(dates, dead2, 'color', 'red', 'd');
legend({"Infected", "Recovered", "Dead"}, 'location', 'northwest');