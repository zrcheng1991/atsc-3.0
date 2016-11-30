%* Revised Sum-Product Algorithm
%*
%* References:
%*   [1] S.J. Johnson, "Low-Density Parity-Check Codes: Design and Decoding", 
%*       Wiley Encyclopedia of Telecommunications, Wiley, Apr. 2003
%* 
%* Author: T.J. Cheng, 2016
%* 
%*   2016-11-30: It works fine for Example 2.6 in [1].
%*

clc
clear all;

TEST = 0;
ITER = 10;

if (~TEST)
    % params
    N       = 16200;
    rate    = 7/15;
    q1      = 24;

    % generate information bits
    s = randi([0 1], 1, N * rate);
    %disp(s);

    % encode information bits
    r = ldpc_enc_B(s, N, rate, q1);

    % generate parity-check matrix
    disp('generating parity-check matrix...');
    H = ldpc_pcmg_B(N, rate, q1);
    fprintf('\tdone\n');
    %spy(H);
else
    r = [ -0.5, 2.5 -4.0 5.0 -3.5 2.5 ];
    H = [ 1 1 0 1 0 0 
          0 1 1 0 1 0 
          1 0 0 0 1 1 
          0 0 1 1 0 1 ];
    H = sparse(H);
end

% sum-product algorithm
fprintf('decoding inner codes...\n');

m = size(H);
M = sparse(m);
E = sparse(m);
L = sparse(1, m(1));

% initialize
for j = 1 : m(1)
    for i = 1 : length(H(j, :))
        if (H(j, i) == 1)
            M(j, i) = r(i);
        end
    end
end

% iteratively decode
for I = 1 : ITER
    for j = 1 : m(1)
        for i = 1 : length(H(j, :))
            if (H(j, i) == 1)
                E(j, i) = 1;
                for ii = 1 : length(H(j, :))
                    if (H(j, ii) == 1 && ii ~= i)
                        E(j, i) = E(j, i) * tanh(M(j, ii) / 2);
                    end
                end
                E(j, i) = 2 * atanh(E(j, i));
            end
        end
    end
    
    for i = 1 : m(2)
        L(i) = r(i);
        for j = 1 : length(H(:, i))
            if (H(j, i) == 1)
                L(i) = L(i) + E(j, i);
            end
        end
        if (L(i) <= 0)
            z(i) = 1;
        else
            z(i) = 0;
        end
    end
    
    if I == ITER || sum(rem(H * z', 2)) == 0
        fprintf('\tdone within %d iterations\n', I);
        break;
    else
        for i = 1 : m(2)
            for j = 1 : length(H(:, i))
                if (H(j, i) == 1)
                    M(j, i) = r(i);
                    for jj = 1 : length(H(:, i))
                        if (H(jj, i) == 1 && jj ~= j)
                            M(j, i) = M(j, i) + E(jj, i);
                        end
                    end
                end
            end
        end
    end
end