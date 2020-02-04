function berout = modem_multipath_ofdm(modemtype, M, N, varargin)
    % MODEM_AWGN Do the modem stuff
    % varargins fs, sps, pathdelays, pathpwr, nsubs, guardbands, prefixlen

    nbits = log2(M);

    % Parse input
    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'fs', 100);
    addParameter(p, 'pathdelays', [0, 2e-9]);
    addParameter(p, 'pathpwr', [0, -10]);
    addParameter(p, 'nsubs', 64);
    addParameter(p, 'guardbands', [6; 6]);
    addParameter(p, 'prefixlen', 16);

    parse(p, varargin{:});

    fs = p.Results.fs;
    pathdelays = p.Results.pathdelays;
    pathpwr = p.Results.pathpwr;
    nsubs = p.Results.nsubs;
    guardbands = p.Results.guardbands;
    prefixlen = p.Results.prefixlen;

    usedsubs = nsubs - length(guardbands);

    nSym = fix(N / nbits / usedsubs);

    % Generate input data
    dataTx = randi([0 M-1], usedsubs, nSym);

    % Modulate the data
    switch modemtype
        case "psk"
            dataMod = pskmod(dataTx, M, 0, "gray");
        case "qam"
            dataMod = qammod(dataTx, M, "gray");
    end % modemtype

    % Do the OFDM modulation
    txSig = ofdmmod(dataMod, nsubs, prefixlen, guardbands);

    % Apply the multipath
    rxSig = apply_multipath(txSig, pathdelays, pathpwr, fs);

    % OFDM Demodulation
    dataRx = ofdmdemod(rxSig, nsubs, prefixlen, prefixlen, guardbands);

    % Finally reverse the modulation
    switch modemtype
        case "psk"
            dataOut = pskdemod(dataRx, M, 0, "gray");
        case "qam"
            dataOut = qamdemod(dataRx, M, "gray");
    end % modemtype

    % Flatten the arrays
    dataTxFlat = reshape(dataTx.', 1, []);
    dataRxFlat = reshape(dataOut.', 1, []);
    rxBits = de2bi(dataRxFlat, nbits);
    txBits = de2bi(dataTxFlat, nbits);
    rxBits = reshape(rxBits.', 1, []);
    txBits = reshape(txBits.', 1, []);

    [~, berout] = biterr(txBits, rxBits);

end % function
