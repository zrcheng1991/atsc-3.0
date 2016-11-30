%* Permutation Matrix
%*
%* References:
%*   [1] K.J. Kim et al., "Low-Density Parity-Check Codes for ATSC 3.0",
%*       IEEE Transactions on Broadcasting, Vol. 59, No. 1, Mar. 2016
%*            
%* Author: T.J. Cheng, 2016
%* 
%*   2016-11-20: It works fine.
%*

function Q = perm_gen_Q(M, L, q)

Q = sparse(M);
for a = 0 : q - 1               % a = i1, j0
    for b = 0 : L - 1           % b = i0, j1
        Q(a * L + b + 1, b * q + a + 1) = 1;
    end
end

return;

