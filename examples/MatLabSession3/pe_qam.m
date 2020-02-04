clear all
clc
M = 16;
Es_No_db = 16;
Pe_symbol = 4*(1-1/sqrt(M))*qfunc(sqrt(3/(M-1)*10^(Es_No_db/10)))
%Assuming Gray code
Pe_bit = Pe_symbol/log2(M)