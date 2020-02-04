% DCS_Assignment_1920
% Rui Oliveira

% ASK FSK PSK QAM

close all;
clear all;

% Parameters
M = 4; % Modulatior order
N = 2^18; % Number of tx bits
tsnr = 0; % Target SNR
clen = 16;
mlen = 2;

% Generate tx stream
bitStream = randi([0 1], 1, N); % Generate data

% Do the FEC encoding
% First we need to adjust the size of the vector
% so that we have one that fit the encoding process
bitStreamFixed = bitStream;

while mod(length(bitStreamFixed), mlen)
    bitStreamFixed = [bitStreamFixed randi([0 1], 1, 1)]; % Fill in a bit
end

% Now the message can fit into the code
% Create a generator polynomial for a cyclic code.
% Create a parity-check matrix by using the generator polynomial.
gpol = cyclpoly(clen, mlen);
parmat = cyclgen(clen, gpol);
% Create a syndrome decoding table by using the parity-check matrix.
trt = syndtable(parmat);
% Encode the data by using the generator polynomial.
encData = encode(bitStreamFixed, clen, mlen, 'cyclic/binary', gpol);

% Now we must fit the message to the PSK symbols
encDataFixed = encData;

while mod(length(encDataFixed), log2(M))
    encDataFixed = [encDataFixed randi([0 1], 1, 1)];
end

% Now that the data fits, we transmit it
databin = reshape(encDataFixed, [length(encDataFixed) / log2(M), log2(M)]);
% Get symbols from the binary stream
data = bi2de(databin);

% Modulate data
if M == 2
    % For M = 2, initial phase is more convenient zero
    txSig = pskmod(data, M, 0, 'gray');
else
    % Else, as per documentation
    txSig = pskmod(data, M, pi / M, 'gray');
end

% Appy AWGN
rxSig = awgn(txSig, tsnr, 'measured');

% Demodulate data
if M == 2
    % For M = 2, initial phase is more convenient zero
    dataOut = pskdemod(rxSig, M, 0, 'gray');
else
    % Else, as per documentation
    dataOut = pskdemod(rxSig, M, pi / M, 'gray');
end

% Now for the recovery process
% First, convert data back into binary
dataRxBin = de2bi(dataOut, log2(M));
dataRxBin = reshape(dataRxBin, 1, []); % This flattens out the array
% Second, adjust the size to the message length (aka, remove what we added)
dataRxBin = dataRxBin(1:length(encData));
% Decode the corrupted sequence
dataRecBinary = decode(dataRxBin, clen, mlen, 'cyclic/binary', gpol, trt);
% Cut the recovered binary into the original shape
dataRecBinary = dataRecBinary(1:length(bitStream));

% Now check what we got
[~, ber] = biterr(bitStream, dataRecBinary);
[~, ser] = symerr(bitStream, dataRecBinary);

% Compare with the expected BER with no FEC
[berteo, serteo] = berawgn(tsnr - 10 * log10(2), "psk", M, 'nondiff');
