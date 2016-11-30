%* Type B (IRA-Structured) LDPC Encoder
%*
%* References:
%*   [1] K.J. Kim et al., "Low-Density Parity-Check Codes for ATSC 3.0",
%*       IEEE Transactions on Broadcasting, Vol. 59, No. 1, Mar. 2016
%*   [2] S.J. Huang., "Implementation and Efficiency Improvment of DVB-T2
%*       Error-Correcting Code", Jun. 2012
%*   [3] Physical Layer Protocol, document ATSC A/322, pp.38-39, ATSC, Sep.
%*       2016     
%* 
%* Author: T.J. Cheng, 2016
%* 
%*   2016-11-20: Error: m exceeds matrix p's dimensions
%*   2016-11-21: Bug fixed (adopts procedure in [2]). It works fine.
%*   2016-11-26: Added parity intlv. described in [3]. It works fine.
%*

function u = ldpc_enc_B(s, N, rate, q1)

% params
K = N * rate;
M = N - K;
L = M / q1;

% load LDPC code matrix
b = ldpc_db(N, rate);

% initialize parity bits
p = zeros(1, M);

% accumulate
for i = 0 : K - 1
    j = floor(i / L);
    z = mod(i, L);
    for l = 0 : sum(b(j + 1, :) >= 0) - 1
        m = mod(b(j + 1, l + 1) + z * q1, M);
        %m = b(j + 1, l + 1) + z;                % here is the bug
        p(m + 1) = xor(p(m + 1), s(i + 1));
    end
end

% sequentially accumulate
for i = 1 : M - 1
    p(i + 1) = xor(p(i + 1), p(i));
end

% concatenate
c = horzcat(s, p);

% interleave parity bits
u = zeros(1, N);
u(:, 1 : K) = c(:, 1 : K);

for s = 0 : L - 1
    for t = 0 : q1 - 1
        u(K + t * L + s + 1) = c(K + s * q1 + t + 1);
    end
end

return;