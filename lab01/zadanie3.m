clc;
clear all;
load('adsl_x.mat');

pref = x(32:2049); 

x = xcorr(x);
x= x/max(x)
plot(x)