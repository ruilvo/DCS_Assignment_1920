function [ber, ser] = does_all_awgn(bitStream, M, snrdb, modtype, fectype, clen, mlen)
    % Do the function for AWGN channels

    % First see if there is any precoding to do
    switch fectype
        case "none"
            bitStreamNew = bitStream;
        case "crc"
            [bitStreamNew, gpol, trt] = do_fec_crc_enc(bitStream, clen, mlen);
    end

    % Reshape the binary stream into the correct size
    % If the bit stream can't be split we can change it to fit
    bitStreamFixed = bitStreamNew;

    while ~(length(bitStreamFixed) / log2(M) == fix(length(bitStreamFixed) / log2(M)))
        bitStreamFixed = [bitStreamFixed randi([0 1], 1, 1)]; % Fill in a bit
    end

    databin = reshape(bitStreamFixed, [length(bitStreamFixed) / log2(M), log2(M)]);
    % Get numbers from that
    data = bi2de(databin);

    switch modtype
        case 'psk'
            % Modulate data
            if M == 2
                % For M = 2, initial phase is more convenient zero
                txSig = pskmod(data, M, 0, 'gray');
            else
                % Else, as per documentation
                txSig = pskmod(data, M, pi / M, 'gray');
            end

        case 'qam'
            txSig = qammod(data, M, 'gray');
    end

    rxSig = awgn(txSig, snrdb, 'measured'); % Appy AWGN

    switch modtype
        case 'psk'
            % Demodulate data
            if M == 2
                % For M = 2, initial phase is more convenient zero
                dataOut = pskdemod(rxSig, M, 0, 'gray');
            else
                % Else, as per documentation
                dataOut = pskdemod(rxSig, M, pi / M, 'gray');
            end

        case 'qam'
            dataOut = qamdemod(data, M, 'gray');
    end

    % Now do the error code recovery
    switch fectype
        case "crc"
            dataRecBinary = do_fec_crc_dec(dataOut, log2(M), clen, mlen, gpol, trt);
            % Now get this binary data again into a proper size
            while mod(length(dataRecBinary), log2(M))
                dataRecBinary(end) = [];
            end

            % And convert back into integer, after reshaping
            dataRecB2D = reshape(dataRecBinary, [length(dataRecBinary) / log2(M), log2(M)]);
            dataOut = bi2de(dataRecB2D);
    end

    % Get parameteres
    csize = min([length(data) length(dataOut)]);
    [~, ber] = biterr(data(1:csize), dataOut(1:csize));
    [~, ser] = symerr(data(1:csize), dataOut(1:csize));
end
