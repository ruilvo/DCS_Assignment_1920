function berout = modem_awgn(modemtype, M, bitsin, tebnos, varargin)
    % MODEM_AWGN Do the modem stuff
    % varargins  fs, sps, rcbeta, rcspan, FECtype, crcCoding, trellisArgs

    N = length(bitsin);
    nbits = log2(M);

    % Parse input
    p = inputParser;
    p.KeepUnmatched = true;
    expectedFECs = {'none', 'cyclic', 'convolutional'};
    addParameter(p, 'fs', 100);
    addParameter(p, 'sps', 16);
    addParameter(p, 'rcbeta', 0.2);
    addParameter(p, 'rcspan', 8);
    addParameter(p, 'FECtype', "none", @(x) any(validatestring(x, expectedFECs)));
    addParameter(p, 'crcCoding', [8, 4]);
    addParameter(p, 'trellisArgs', {[5 4], [23 35 0; 0 5 13]});

    parse(p, varargin{:});

    fs = p.Results.fs;
    sps = p.Results.sps;
    rcbeta = p.Results.rcbeta;
    rcspan = p.Results.rcspan;
    FECtype = p.Results.FECtype;
    crcCoding = p.Results.crcCoding;
    trellisArgs = p.Results.trellisArgs;

    switch FECtype
        case "none"
            bitStream = bitsin;
        case "cyclic"
            clen = crcCoding(1);
            mlen = crcCoding(2);

            % Do the FEC encoding
            % First we need to adjust the bitstream so that it fits
            % We need to add mlen-mod(len,mlen) bits
            bitStreamFEC = [bitsin; randi([0 1], mlen - mod(N, mlen), 1)];

            % Now the message can fit into the code
            % Create a generator polynomial for a cyclic code.
            gpol = cyclpoly(clen, mlen);
            % Create a parity-check matrix by using the generator polynomial.
            parmat = cyclgen(clen, gpol);
            % Create a syndrome decoding table by using the parity-check matrix.
            trt = syndtable(parmat);

            % Encode the data by using the generator polynomial.
            bitStream = encode(bitStreamFEC, clen, mlen, 'cyclic/binary', gpol);
        case "convolutional"
            % Convert convolutional code polynomials to trellis description
            switch length(trellisArgs)
                case 2
                    trellis = poly2trellis(trellisArgs{1}, trellisArgs{2});
                case 3
                    trellis = poly2trellis(trellisArgs{1}, trellisArgs{2}, trellisArgs{3});
            end % trellisArgs

            % Convolutionally encode the data
            bitStream = convenc(bitsin, trellis);

            bsupsampling = log2(trellis.numOutputSymbols); % /trellis.numInputSymbols;
        otherwise
            bitStream = bitsin;
    end % FECtype

    % First, fit the array into the symbols
    % We need to add nbits-mod(len,nbits) bits
    symbstuffn = nbits - mod(length(bitStream), nbits);
    bitStreamFixed = [bitStream; randi([0 1], symbstuffn, 1)];

    % We reshape the stream to convert into symbols
    databin = reshape(bitStreamFixed, [nbits, length(bitStreamFixed) / nbits]);
    databin = databin.';
    % Get symbols from the binary stream
    data = bi2de(databin);

    % Modulate the data
    switch modemtype
        case "psk"

            if M == 2
                % For M = 2, initial phase is more convenient zero
                txSigUf = pskmod(data, M, 0, 'gray');
            else
                % Else, as per documentation
                txSigUf = pskmod(data, M, pi / M, 'gray');
            end

            usefilter = true;
        case "qam"
            txSigUf = qammod(data, M, 'gray');
            usefilter = true;
        case "fsk"
            fsep = fs / 2 / M;
            txSig = fskmod(data, M, fsep, sps, fs, "cont", 'gray');
            usefilter = false;
    end % modemtype

    % Apply raised cossine filter is need be
    if usefilter
        % Design matched filters
        txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor', rcbeta, ...
            'FilterSpanInSymbols', rcspan, 'OutputSamplesPerSymbol', sps);
        rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor', rcbeta, ...
            'FilterSpanInSymbols', rcspan, 'InputSamplesPerSymbol', sps, ...
            'DecimationFactor', sps);

        % Filter data on tx
        txSig = txfilter(txSigUf);
    end % usefilter

    % Now we are ready to iterate over the various SNRs
    berout = zeros(length(tebnos), 1);

    for a = 1:length(tebnos)
        tsnr = tebnos(a) + 10 * log10(nbits) - 10 * log10(sps);

        fomatSpec = "SNR = %.3e; %s, M = %d\n";
        fprintf(fomatSpec, tsnr, modemtype, M);
        drawnow;

        % Appy AWGN
        rxSig = add_awgn_noise(txSig, tsnr);

        % If needed, apply the reverse filter
        if usefilter
            % If the filter was applied, apply the matched filter (and downsampling)
            rxSigUf = rxfilter(rxSig);
            rxSigUf = rxSigUf(rcspan + 1:end);
        end % usefilter

        % Now demodulate
        switch modemtype
            case "psk"

                if M == 2
                    % For M = 2, initial phase is more convenient zero
                    dataOut = pskdemod(rxSigUf, M, 0, 'gray');
                else
                    % Else, as per documentation
                    dataOut = pskdemod(rxSigUf, M, pi / M, 'gray');
                end

            case "qam"
                dataOut = qamdemod(rxSigUf, M, 'gray');
            case "fsk"
                dataOut = fskdemod(rxSig, M, fsep, sps, fs, "gray");
        end % modemtype

        % Convert symbols to bits
        rxbits = de2bi(dataOut, nbits);
        % Convert to line array
        rxbits = reshape(rxbits.', 1, []);
        rxbits = rxbits.';

        % Now we must be careful with decoding
        if FECtype == "none"
            % Nothing special to do.
            % Just make sure arrays have the same length
            % This automatically deals with not having
            % the last couple bits in the usefilter case
            cutlen = min(length(rxbits), length(bitsin));

            [~, ber] = biterr(rxbits(1:cutlen), bitsin(1:cutlen));

        elseif FECtype == "cyclic"

            % First, trim the data into a length that fits the mlen
            trimsize = mod(length(rxbits), clen);

            dataRxBin = rxbits(1:end - trimsize);

            % Decode the corrupted sequence
            dataRecBinary = decode(dataRxBin, clen, mlen, 'cyclic/binary', gpol, trt);

            % Cut the array into it's maximum useful length
            cutlen = min(length(dataRecBinary), length(bitsin));

            [~, ber] = biterr(dataRecBinary(1:cutlen), bitsin(1:cutlen));

        elseif FECtype == "convolutional"
            % First, trim the data into a length that fits the codes
            trimsize = mod(length(rxbits), bsupsampling);

            dataRxBin = rxbits(1:end - trimsize);

            % Decode the data
            dataRecBinary = vitdec(dataRxBin, trellis, 34, 'trunc', 'hard');

            % Trim to only data that makes sense
            cutlen = min(length(dataRecBinary), length(bitsin));

            [~, ber] = biterr(dataRecBinary(1:cutlen), bitsin(1:cutlen));
        end % FECtype

        berout(a) = ber;

    end % tebnos

end % function
