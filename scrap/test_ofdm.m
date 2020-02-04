% Testing ODFM (de)modulation
% With PSK symbols

close all;
clear all;

%% Parameters
% Number of tx bits
N = 2^15;

% Modulation order
M = 4;

% OFDM parameters
nsubs = 64;
guardbands = [6; 6];
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
usedsubs = nsubs - sum(guardbands);
frameSize = nbits * usedsubs;

% The bitsin must fit within K frames
Neff = N + frameSize - mod(N, frameSize);

% Generate input data
bitsin = randi([0 1], Neff, 1);

% Create modems

% Create the ODFM modem
ofdmMod = comm.OFDMModulator("FFTLength", nsubs, "NumGuardBandCarriers", guardbands, ...
    "CyclicPrefixLength", prefixlen);
ofdmDemod = comm.OFDMDemodulator("FFTLength", nsubs, "NumGuardBandCarriers", guardbands, ...
    "CyclicPrefixLength", 16);
ofdmDims = info(ofdmMod);

% Create the PSK modem
pskMod = comm.PSKModulator("ModulationOrder", M, "PhaseOffset", 0, ...
    "BitInput", true, "SymbolMapping", "Gray");
pskDemod = comm.PSKDemodulator("ModulationOrder", M, "PhaseOffset", 0, ...
    "BitOutput", true, "SymbolMapping", "Gray", "DecisionMethod", "Hard decision");

% Apply PSK modulation
dataTx = pskMod(bitsin);

% Reshape symbols into frames
nframes = length(dataTx) / usedsubs;
dataTxFramed = reshape(dataTx, [usedsubs, nframes]);

% Apply OFDM modulation;
% We must do this frame by frame
txSig = [];

for k = 1:nframes
    tempSlice = ofdmMod(dataTxFramed(:, k));
    txSig = [txSig; tempSlice];
end

% TODO: do something with the channel (like delays and/or AWGN)
rxSig = txSig;

% Now to demodulate,
% first, reshape rxSig to an approriate shape
rxSigFramed = reshape(rxSig, ofdmDims.OutputSize(1), nframes);

% Apply OFDM demod
% We must do this frame by frame
dataRx = [];

for k = 1:nframes
    tempSlice = ofdmDemod(rxSigFramed(:, k));
    dataRx = [dataRx; tempSlice]; % Note the ;
end

% Finally apply PSK demod
dataOut = pskDemod(dataRx);

[number, ratio] = biterr(bitsin, dataOut);

scope = dsp.SpectrumAnalyzer('SampleRate',20e6);
scope(rxSig);
