clear all;
close all;

% Number of tx bits
N = 128;
% Samples per symbol
sps = 32;
% Sample frequency
fs = 1e6; % Hz
% RC filter parameters
rcbeta = 0.2;
rcspan = 10;

tsnr =100;

M = 4;

% Generate tx stream
bitStream = randi([0 1], N, 1);

nbits = log2(M);

% First, fit the array into the symbols
% We need to add nbits-mod(len,nbits) bits
bitStreamFixed = [bitStream; randi([0 1], nbits - mod(N, nbits), 1)];

% We reshape the stream to convert into symbols
databin = reshape(bitStreamFixed, [length(bitStreamFixed) / nbits, nbits]);
% Get symbols from the binary stream
data = bi2de(databin);

txSigUf = pskmod(data, M, 0, 'gray');

txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor', rcbeta, ...
    'FilterSpanInSymbols', rcspan, 'OutputSamplesPerSymbol', sps);

rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor', rcbeta, ...
    'FilterSpanInSymbols', rcspan, 'InputSamplesPerSymbol', sps, ...
    'DecimationFactor', sps);

txSig = txfilter(txSigUf);

rxSig = awgn(txSig, tsnr, 'measured');

rxSigUf = rxfilter(rxSig);

dataOut = pskdemod(rxSigUf(rcspan + 1:end), M, 0, 'gray');

dataRxBin = de2bi(dataOut, nbits);
dataRxBin = reshape(dataRxBin, [], 1);

datacut = data(1:end - rcspan);
datacutBin = de2bi(datacut, nbits);
datacutBin = reshape(datacutBin, [], 1);

figure(1); hold on;
plot(real(dataOut));
plot(real(datacut));

figure(2); hold on;
plot(real(dataRxBin));
plot(real(datacutBin));


