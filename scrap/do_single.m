% DCS class assigment, MAP-tele, 2019/2020
% Rui Oliveira

close all;
clear all;


% Parameters
M = 4;
N = 2^18; % Number of tx bits. Must be divisible by 6 (because log2(64))

tsnr = 0;

modtype = "psk";
fectype = "crc";
clen = 16;
mlen = 2;


% Generate tx stream
bitStream = randi([0 1],1,N); % Generate data

[ber,ser] = does_all_awgn(bitStream,M,tsnr,modtype,fectype,clen,mlen);

switch modtype
    case 'psk'
        [berteo, serteo] = berawgn(tsnr - 10*log10(2),modtype,M, 'nondiff'); 
    case 'qam'
        [berteo, serteo] = berawgn(tsnr - 10*log10(2),modtype,M); 
end
