%* Type A (MET-Structured) LDPC Encoder
%*
%* References:
%*   [1] K.J. Kim et al., "Low-Density Parity-Check Codes for ATSC 3.0",
%*       IEEE Transactions on Broadcasting, Vol. 59, No. 1, Mar. 2016
%*   [2] S.I. Park et al., US Patent No. 20160049954 A1, 2016
%* 
%* Author: T.J. Cheng, 2016
%*
%*   2016-11-24: checked by Dr. You
%*

function c = ldpc_enc_A(s, N, rate, M1, M2, q1, q2)

% params
K = N * rate;
L = M1 / q1;

% load LDPC code matrix
b = ldpc_db(N, rate);

% 1. initialize
c = zeros(1, N);
c(1 : K) = s;
p = zeros(1, M1 + M2);

% 2. accumulate
for i = 0 : K - 1
    j = floor(i / L);
    for l = 0 : sum(b(j + 1, :) >= 0) - 1
        x = b(j + 1, l + 1);
        if x < M1
            x = mod(x + i * q1, M1);
        else
            x = M1 + mod(x - M1 + i * q2, M2);
        end
        p(x + 1) = xor(p(x + 1), s(i + 1));
    end
end

% 3. sequentially accumulate
for i = 2 : M1
    p(i) = xor(p(i), p(i - 1));
end

% 4. interleave p1
for t = 0 : q1 - 1
    for s = 0 : L - 1
        c(K + t * L + s + 1) = p(s * q1 + t + 1);
    end
end

% 5. accumulate p1 to p2
for i = K : K + M1 - 1
    j = floor(i / L);
    for l = 0 : sum(b(j + 1, :) >= 0) - 1
        x = b(j + 1, l + 1);
        if x < M1
            x = mod(x + i * q1, M1);
        else
            x = M1 + mod(x - M1 + i * q2, M2);
        end
        p(x + 1) = xor(p(x + 1), c(i + 1));
    end
end

% for i = 2 : M1 + M2
%     p(i) = xor(p(i), p(i - 1));
% end
 
% 6. interleave p2
for t = 0 : q2 - 1
    for s = 0 : L - 1
        c(K + M1 + t * L + s + 1) = p(M1 + s * q2 + t + 1);
    end
end

return;