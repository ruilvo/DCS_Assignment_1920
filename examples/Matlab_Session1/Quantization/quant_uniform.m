clear all
close all
clc
bdclose('all');

R=1000;
fs=10*R;
N = 8
L = 2^N;
m_max = 1;
m_min = -m_max;
delta = (m_max - m_min)/L;
boundary = m_min : delta : m_max;
codebook = m_min+delta/2 : delta : m_max-delta/2;
atenuation_dB = 20
% 4.902 normalizes the amplitude of the sum of thtwo sine waves
k = 1/4.902/(sqrt(10^(atenuation_dB/10)));
open_system('quantizer_uniform');
sim('quantizer_uniform');

% theoretical value uniform quantizer
NQ_t = delta^2/12;

Amp1 = 3;
Amp2 = 2;
S_t = (Amp1*k)^2/2 + (Amp2*k)^2/2;
S_NQ_t = S_t/NQ_t;
S_NQ_t_db = 10*log10(S_NQ_t)

% experimental values
NQ = mean(eq.^2);
S = mean(x.^2);
S_NQ = S/NQ;
S_NQdB = 10*log10(S_NQ)