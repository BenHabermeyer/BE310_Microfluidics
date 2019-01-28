%BE310 STAT Lecture

%% error 
clc
T = [65.1, 64.8, 65.0, 65.1, 64.9, 65.1];
estimate = mean(T) %estimate is mean
sd = std(T); %uncertainty is standard deviation
uncertainty = sd / sqrt(6)

%follow up
E_kb = 1 * ones(size(T)); %Kelvin
A = 1;
%convert to kelvin
T = T + 273.15;

k = A .* exp(E_kb ./ T);
estimate_k = mean(k)
%take the derivative and multiply by uncertainty in T
uncertainty_k = A * E_kb(1)  * estimate^(-2) * exp(-E_kb(1) / estimate) * uncertainty

%% flow cytometer
clc
close all
data = table2array(flowcytometerdata);

cellA = data(:,1);
cellB = data(:,2);

figure;
subplot(2,1,1);
histogram(cellA);
title('CellA EpCAM distribution N=1000');
xlabel('count');
subplot(2,1,2);
histogram(cellB);
title('CellB EpCAM distribution N=1000');
ylabel('Fluorescence intensity');
xlabel('count');

%estimates
N = 1000;
estimateA = mean(cellA)
uncertaintyA = std(cellA) / sqrt(N)
estimateB = mean(cellB)
uncertaintyB = std(cellB) / sqrt(N)


estimate_difference = estimateA - estimateB
uncertainty_difference = sqrt(uncertaintyA^2 + uncertaintyB^2)

%there is a difference between the two distributions



