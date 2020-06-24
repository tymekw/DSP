clear all;
close all;

As = 40;
Ap =3;

Fs1 = 96*10^6 - (120*10^3);
Fp1 = 96*10^6 - (100*10^3);

Fs2 = 96*10^6 + (120*10^3);
Fp2 = 96*10^6 + (100*10^3);

D = fdesign.bandpass(Fs1,Fp1,Fp2,Fs2,As,Ap,As,40e7);
H = design(D,'ellip');
fvtool(H);