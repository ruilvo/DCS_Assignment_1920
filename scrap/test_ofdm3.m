% Testing ODFM (de)modulation
% With PSK symbols

close all;
clear all;

%% Parameters
% Number of tx bits
N = 2^15;

% Modulation order
M = 16;

fs = 20e6;

% OFDM parameters
nsubs = 64;
guardbands = [1:6 nsubs-5:nsubs]';
prefixlen = 16;
% Note on the assignment
% N = 64 - 6 - 6 = 52 useful subcariers, as intended
% 0.8us/(4us-0.8us) = 0.25
% 0.25 * 64 = 16, therefore 16 for the prefix length
% We don't have a notion of time here.
% It'll probably be important to find out how to select
% fs such as the symbol timings yield correct values.

%% Simulation
% Computed parameters
nbits = log2(M);
usedsubs = nsubs - length(guardbands);
frameSize = nbits * usedsubs;

nSym = fix(N / nbits / usedsubs);

% Generate input data
dataTx = randi([0 M-1], usedsubs, nSym);

dataMod = qammod(dataTx, M, "gray","UnitAveragePower",true);
txSig = ofdmmod(dataMod, nsubs, prefixlen, guardbands);

% Simulate multipath

pathdelays = [0, 1.2e-6];
pathpwr = [0, -6];

rxSig = apply_multipath(txSig, pathdelays, pathpwr, fs);

dataRx = ofdmdemod(rxSig, nsubs, prefixlen, prefixlen, guardbands);

% Finally apply PSK demod
dataOut = qamdemod(dataRx, M, "gray","UnitAveragePower",true);

% Flatten the arrays
dataTxFlat = reshape(dataTx.', 1, []);
dataRxFlat = reshape(dataOut.', 1, []);
rxBits = de2bi(dataRxFlat, nbits);
txBits = de2bi(dataTxFlat, nbits);
rxBits = reshape(rxBits.', 1, []);
txBits = reshape(txBits.', 1, []);

[number, ratio] = biterr(txBits, rxBits)

%scope = dsp.SpectrumAnalyzer('SampleRate', fs);
%scope(txSig);

