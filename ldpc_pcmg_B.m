%* Type B (IRA-Structured) LDPC Parity-Check Matrix Generator
%*
%* References:
%*   [1] K.J. Kim et al., "Low-Density Parity-Check Codes for ATSC 3.0",
%*       IEEE Transactions on Broadcasting, Vol. 59, No. 1, Mar. 2016
%*   [2] Physical Layer Protocol, document ATSC A/322, pp.38-39, ATSC, Sep.
%*       2016 
%*            
%* Author:  T.J. Cheng, 2016
%* 
%*   2016-11-20: It works fine (without converting to QC-form).
%*   2016-11-26: Added parity intlv. described in [2]. It works fine.
%*

function H = ldpc_pcmg_B(N, rate, q1)

% params
K = N * rate;
M = N - K;
L = M / q1;
k = K / L;

% load LDPC code matrix
b = ldpc_db(N, rate);

% build information part
Hi = sparse(M, K);
for j = 0 : k - 1
    for l = 0 : sum(b(j + 1, :) >= 0) - 1
        for z = 0 : L - 1
            Hi(mod(b(j + 1, l + 1) + z * q1, M) + 1, j * L + z + 1) = 1;
        end
    end
end

% build parity part
Hp = speye(M);
Hp = spdiags(ones(M - 1, 1), -1, Hp);

% concatenate two parts
H = horzcat(Hi, Hp);

% converting to QC-form
Q = perm_gen_Q(M, L, q1);
P = zeros(N);
P(1 : K, 1 : K) = speye(K);
P(K + 1 : end, K + 1 : end) = Q';

H = Q * H * P;

return;