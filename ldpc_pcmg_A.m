%* Type A (MET-Structured) LDPC Parity-Check Matrix Generator
%*
%* References:
%*   [1] K.J. Kim et al., "Low-Density Parity-Check Codes for ATSC 3.0",
%*       IEEE Transactions on Broadcasting, Vol. 59, No. 1, Mar. 2016
%*            
%* Author:  T.J. Cheng, 2016
%*
%*   2016-11-20: It works fine.
%*

function H = ldpc_pcmg_A(N, rate, M1, M2, q1, q2)

% params
K = N * rate;
L = M1 / q1;

% load LDPC code matrix
b = ldpc_db(N, rate);

% build information part (matrix A)
H = sparse(N - K, N);
k = K / L;
for j = 0 : k - 1
    for l = 0 : sum(b(j + 1, :) >= 0) - 1
        for z = 0 : L - 1
            x = b(j + 1, l + 1);
            if x < M1
                H(mod(b(j + 1, l + 1) + z * q1, M1) + 1, j * L + z + 1) = 1;
            end
        end
    end
end

% build information and parity-1 part (matrix C and D')
k = (K + M1) / L;
for j = 0 : k - 1
    for l = 0 : sum(b(j + 1, :) >= 0) - 1
        for z = 0 : L - 1
            x = b(j + 1, l + 1);
            if x >= M1
                H(M1 + mod(b(j + 1, l + 1) - M1 + z * q2, M2) + 1, j * L + z + 1) = 1;
            end
        end
    end
end

% build matrix B (dual diagonal matrix)
Hb = speye(M1);
Hb = spdiags(ones(M1 - 1, 1), -1, Hb);
H(1 : M1, K + 1 : K + M1) = Hb;

% build matrix E (identity matrix)
He = speye(M2);
H(M1 + 1 : M1 + M2, K + M1 + 1 : N) = He;

% permutate matrix A and B
Q1 = perm_gen_Q(M1, L, q1);
Ha = H(1 : M1, 1 : K);
H(1 : M1, 1 : K) = Q1 * Ha;
Hb = H(1 : M1, K + 1 : K + M1);
H(1 : M1, K + 1 : K + M1) = Q1 * Hb * Q1';

% permutate matrix C and D'
Q2 = perm_gen_Q(M2, L, q2);
Hc = H(M1 + 1 : M1 + M2, 1 : K);
H(M1 + 1 : M1 + M2, 1 : K) = Q2 * Hc;
Hd = H(M1 + 1 : M1 + M2, K + 1 : K + M1);
H(M1 + 1 : M1 + M2, K + 1 : K + M1) = Q2 * Hd;

return;