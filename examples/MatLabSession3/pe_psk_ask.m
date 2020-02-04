clear all
close all
clc

%bpsk theoretical pe
A=1;
Tb=1/1000;
Eb=A^2*Tb/2;
var=10;
fs=100/Tb;
No=2*var/fs;
Pe_psk = qfunc(sqrt(2*Eb/No))


%ask theoretical pe
A=1;
Tb=1/1000;
Eb=A^2*Tb/2;
var=2.5;
No=2*var/fs;
Pe_ask = qfunc(sqrt(Eb/No/2))
