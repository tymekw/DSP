clc;
clear all;

[data,fs]=audioread('lab0cps.m4a');
%fs czestotliwosc probkowania
%sound(data, fs);

% 1 sluchawka
samples1 = data(:,1);
sample_len1 = length(samples1)/fs;
sredenia1 = mean(samples1);
odchylenie1 = std(samples1);
max1 = max(samples1);
min1 = min(samples1);

newFs=3000;
%sound(data,newFs);


%2 sluchawka
samples2 = data(:,2);
sample_len2 = length(samples2)/fs;
sredenia2 = mean(samples2);
odchylenie2 = std(samples2);
max2 = max(samples2);
min2 = min(samples2);

sound(data,30000);

plot(samples1,'b');
hold on;
plot(samples2,'r');