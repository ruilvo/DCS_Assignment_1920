clear all
close all
clc
bdclose('all');


[x,fs]=audioread('speech.wav');
info = audioinfo('speech.wav');
%bits_per_sample = info.BitsPerSample; 
N = 8
L = 2^N;
m_max = 1;
m_min = -m_max;
delta = (m_max - m_min)/L;
boundary = m_min : delta : m_max;
codebook = m_min+delta/2 : delta : m_max-delta/2;
atenuation_dB = 0
k=sqrt(10^(atenuation_dB/10));
A=87.6;

open_system('speech_uniform');
sim('speech_uniform');


NQ = mean(eq.^2);
S = mean(xout.^2);
S_NQ = S/NQ;


S_NQdB = 10*log10(S_NQ)

sound(k*xq,fs);

