function berout = modem_multipath(modemtype, M, bitsin, varargin)
    % MODEM_AWGN Do the modem stuff
    % varargins  fs, sps, rcbeta, rcspan, fsep, pathdelays, pathpwr

    N = length(bitsin);
    nbits = log2(M);

    % Parse input
    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'fs', 100);
    addParameter(p, 'sps', 16);
    addParameter(p, 'rcbeta', 0.2);
    addParameter(p, 'rcspan', 8);
    addParameter(p, 'fsep', 10);
    addParameter(p, 'pathdelays', [0, 2e-9]);
    addParameter(p, 'pathpwr', [0, -10]);

    parse(p, varargin{:});

    fs = p.Results.fs;
    sps = p.Results.sps;
    rcbeta = p.Results.rcbeta;
    rcspan = p.Results.rcspan;
    fsep = p.Results.fsep;
    pathdelays = p.Results.pathdelays;
    pathpwr = p.Results.pathpwr;

    % First, fit the array into the symbols
    % We need to add nbits-mod(len,nbits) bits
    symbstuffn = nbits - mod(length(bitsin), nbits);
    bitStreamFixed = [bitsin; randi([0 1], symbstuffn, 1)];

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

    % Appy multipath
    rxSig = apply_multipath(txSig, pathdelays, pathpwr, fs);

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
    % Just make sure arrays have the same length
    % This automatically deals with not having
    % the last couple bits in the usefilter case
    cutlen = min(length(rxbits), length(bitsin));

    [~, berout] = biterr(rxbits(1:cutlen), bitsin(1:cutlen));

end % function
