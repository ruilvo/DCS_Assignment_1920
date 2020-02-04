function [encData, gpol, trt] = do_fec_crc_enc(bitStream, n, k)
    % Transforms an array of bits into a one with error correction

    % As with before, we'll fit the bitstream to the correct size
    % It must be splittable by k
    bitStreamFixed = bitStream;

    while ~(mod(length(bitStreamFixed), k))
        bitStreamFixed = [bitStreamFixed randi([0 1], 1, 1)]; % Fill in a bit
    end

    % Now the message can fit into the code
    % Create a generator polynomial for a cyclic code.
    % Create a parity-check matrix by using the generator polynomial.
    gpol = cyclpoly(n, k);
    parmat = cyclgen(n, gpol);
    % Create a syndrome decoding tabled by using the parity-check matrix.
    trt = syndtable(parmat);
    % Encode the data by using the generator polynomial.
    encData = encode(bitStreamFixed, n, k, 'cyclic/binary', gpol);
end
