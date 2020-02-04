close all
clear all
clc

%Detection at the centre Theoretical - AMI Case
A=1;
threshold=A/2;
var_1=0.05;
Pe_centre_ami=1.5*qfunc(threshold/sqrt(var_1)) - 0.5*(qfunc(3*threshold/sqrt(var_1)))


%Detection with Match Filter Theoretical - AMI Case
Tb=1e-3;
A=1;
Eb=A^2*Tb;
fs=100/Tb;
var_2=5;
No=2*var_2/fs;
Pe_match_ami=1.5*qfunc(sqrt(2*Eb/4/No)) - .5*qfunc(sqrt(2*9*Eb/4/No))