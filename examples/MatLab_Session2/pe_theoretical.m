close all
clear all
clc

%Detection at the centre Theoretical - Polar Case
A=1;
delta_v=2*A;
var_1=0.2
Pe_centre = qfunc(delta_v/(2*sqrt(var_1)))
%Pe_centre = 0.5*erfc((delta_v/(2*sqrt(2*var_1))))

%Matched Filter Theoretical - Polar Case
Tb=1e-3;
A=1;
Eb=A^2*Tb;
fs=100/Tb;
var_2=20
No=2*var_2/fs;
Pe_match = qfunc(sqrt(2*Eb/No))
%Pe_match = 0.5*erfc(sqrt(Eb/No))

