% DCS class assigment, MAP-tele, 2019/2020
% Rui Oliveira

close all;
clear all;


% Parameters
Ms = [2 4 16 64];
N = 2^17*3; % Number of tx bits. Must be divisible by 6 (because log2(64))

% Parameters target
minsnr = -20;
maxsnr = 20;
snrstep = 2;

tsnrs = minsnr:snrstep:maxsnr;

Mmodtypes = ["psk", "qam"];

% Results containers 
resultsber = zeros(length(tsnrs),length(Ms)*length(Mmodtypes));
resultsser = zeros(length(tsnrs),length(Ms)*length(Mmodtypes));
resultsberteo = zeros(length(tsnrs),length(Ms)*length(Mmodtypes));
resultsserteo = zeros(length(tsnrs),length(Ms)*length(Mmodtypes));

% Generate tx stream
bitStream = randi([0 1],1,N); % Generate data


for k = 1:length(tsnrs)
    tsnr = tsnrs(k)
    
    v = 0;
    for modtype = Mmodtypes
        for l = 1:length(Ms)
            M = Ms(l);
            
            % There is no QAM below 4
            if modtype == "qam" & M<4
                break
            end
            [ber,ser] = does_all(bitStream,M,tsnr,modtype);

            resultsber(k,v+l) = ber;
            resultsser(k,v+l) = ser;
            
            switch modtype
                case 'psk'
                    [berteo, serteo] = berawgn(tsnr - 10*log10(2),modtype,M, 'nondiff'); 
                case 'qam'
                    [berteo, serteo] = berawgn(tsnr - 10*log10(2),modtype,M); 
            end
            
            resultsberteo(k,v+l) = berteo;
            resultsserteo(k,v+l) = serteo;
        end
        v = v+1;
    end

end

figure(1);
title("BER");
plot(tsnrs, resultsber, '*'); hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(tsnrs, resultsberteo, '--'); hold off;
figure(2);
title("SER");
plot(tsnrs, resultsser, '*'); hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(tsnrs, resultsserteo, '--'); hold off;
