
clear all; close all; clc;

%% Generate pilot
fp = 19e3;                  % Frequency
fs = 8e3;                   % Sampling frequency
t = 0:1/fs:1-1/fs;          % Time vector
p1 = sin(2*pi*fp*t + pi/4); % Pilot

%% Estimate pilot's frequency and phase
freq = 2*pi*fp/fs;
theta = zeros(1,length(p1)+1);
alpha = 1e-2;
beta = alpha^2/4;
for n = 1 : length(p1)
    perr = -p1(n)*sin(theta(n));
    theta(n+1) = theta(n) + freq + alpha*perr;
    freq = freq + beta*perr;
end

error = 2*pi*fp*t-pi/4 - theta(1:end-1);

figure;
%plot(t, 2*pi*fp*t, 'b', t, theta(1:end-1), 'r'); xlabel('t [s]');
%title('Pilot frequency'); legend('Original const. phase', 'Estimated const. phase');
plot(t, error, 'k'); title('Frequency error - const. phase');

%% Generate pilot with FM
df = 10;    % Modulation depth
fm = 0.1;   % Modulation frequency
m = df*sin(2*pi*fm*t);  % Modulating signal
p2 = sin(2*pi*fp*t+m);  % Pilot with FM

%% Estimate pilot's frequency and phase
freq = 2*pi*fp/fs;
theta = zeros(1,length(p2)+1);
alpha = 1e-2;
beta = alpha^2/4;
for n = 1 : length(p2)
    perr = -p2(n)*sin(theta(n));
    theta(n+1) = theta(n) + freq + alpha*perr;
    freq = freq + beta*perr;
end

error = 2*pi*fp*t+m-pi/2 - theta(1:end-1);

figure;
%plot(t, 2*pi*fp*t+m, 'b', t, theta(1:end-1), 'r'); xlabel('t [s]');
%title('Pilot frequency'); legend('Original FM', 'Estimated FM');
plot(t, error, 'k'); title('Frequency error - FM');

%% Add noise to pilot
%p1n = awgn(p1, 0, 'measured');
%p1n = awgn(p1, 5, 'measured');
%p1n = awgn(p1, 10, 'measured');
p1n = awgn(p1, 20, 'measured');

freq = 2*pi*fp/fs;
theta = zeros(1,length(p1n)+1);
alpha = 1e-2;
beta = alpha^2/4;
for n = 1 : length(p1n)
    perr = -p1n(n)*sin(theta(n));
    theta(n+1) = theta(n) + freq + alpha*perr;
    freq = freq + beta*perr;
end

error = 2*pi*fp*t+pi/4 - theta(1:end-1);

figure;
plot(t, error, 'k'); xlabel('t [s]');
title('Frequency error');
