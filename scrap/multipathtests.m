% DCS_Assignment_1920
% Rui Oliveira

close all;
clear all;

%--- Parameters ---
% Number of tx bits
N = 2^18;
% Samples per symbol
sps = 128;
% Sample frequency
fs = 1e4; % Hz
% RC filter parameters
rcbeta = 0.2;
rcspan = 10;

% PSK parameters
M = 2;

% Generate tx stream
bitsin = randi([0 1], N, 1);
symbstuffn = nbits - mod(length(bitStream), nbits);
bitStreamFixed = [bitsin; randi([0 1], symbstuffn, 1)];

% We reshape the stream to convert into symbols
databin = reshape(bitStreamFixed, [nbits, length(bitStreamFixed) / nbits]);
databin = databin.';

% Get symbols from the binary stream
data = bi2de(databin);

txSigUf = pskmod(data, M, 0, 'gray');

% Design matched filters
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor', rcbeta, ...
'FilterSpanInSymbols', rcspan, 'OutputSamplesPerSymbol', sps);
rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor', rcbeta, ...
'FilterSpanInSymbols', rcspan, 'InputSamplesPerSymbol', sps, ...
'DecimationFactor', sps);

% Filter data on tx
txSig = txfilter(txSigUf);
