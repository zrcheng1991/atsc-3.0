%* Sum-Product Algorithm
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

TEST = 1;
ITER = 5;

if (~TEST)
    % params
    N       = 16200;
    rate    = 7/15;
    q1      = 24;

    % generate information bits
    s = randi([0 1], 1, N * rate);
    %disp(s);

    % encode information bits
    c = ldpc_enc_B(s, N, rate, q1);

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
end

A = [ 1 3; 1 2; 2 4; 1 4; 2 3; 3 4 ];
B = [ 1 2 4; 2 3 5; 1 5 6; 3 4 6 ];

fprintf('decoding inner codes...\n');

% initialize
m = size(H);
M = zeros(m);
for j = 1 : m(1)
    for p = 1 : sum(B(j, :) > 0)
        i = B(j, p);
        M(j, i) = r(i);
    end
end

M

% iteratively decode
for I = 1 : ITER
    for j = 1 : m(1)
        for p = 1 : sum(B(j, :) > 0)
            i = B(j, p);
            E(j, i) = 1;
            for q = 1 : sum(B(j, :) > 0)
                ii = B(j, q);
                if (ii ~= i)
                    E(j, i) = E(j, i) * tanh(M(j, ii) / 2);
                end
            end
            E(j, i) = 2 * atanh(E(j, i));
        end
    end
    
    for i = 1 : m(2)
        L(i) = r(i);
        for p = 1 : sum(A(i, :) > 0)
            j = A(i, p);
            L(i) = L(i) + E(j, i);
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
            for p = 1 : sum(A(i, :) > 0)
                j = A(i, p);
                M(j, i) = r(i);
                for q = 1 : sum(A(i, :) > 0)
                    jj = A(i, q);
                    if (jj ~= j)
                        M(j, i) = M(j, i) + E(jj, i);
                    end
                end
            end
        end
    end
end