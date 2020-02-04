% Testing ODFM (de)modulation
% With PSK symbols

% close all;
% clear all;

%% Parameters
% Number of tx bits
N = 2^10;

% Modulation order
M = 4;

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

timedelays = linspace(0, 4e-6, 400);
scndchanpwr = -3; % dB

testsDelays = {};

for a = 1:length(timedelays)
    testsDelays{a} = {[0 timedelays(a)], [0 scndchanpwr]};
end

fs = 20e6;

results = zeros(length(testsDelays), 1);

for i = 1:length(results)
    
    fprintf("%d / %d\n", i, length(results));

    pathdelays = testsDelays{i}{1};
    pathpwr = testsDelays{i}{2};

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

    [number, ratio] = biterr(txBits, rxBits);

    results(i) = ratio;
    
end % for results

%scope = dsp.SpectrumAnalyzer('SampleRate',20e6);
%scope(txSig);

figure(1);
plot(timedelays, results);
