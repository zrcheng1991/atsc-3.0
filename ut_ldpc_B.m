%* Unit Test for Type B (IRA-Structured) LDPC Encoder
%*
%* References:
%*   [1] K.J. Kim et al., "Low-Density Parity-Check Codes for ATSC 3.0",
%*       IEEE Transactions on Broadcasting, Vol. 59, No. 1, Mar. 2016
%*   [2] MathWorks, "Error Detection and Correction", Available: 
%*         https://www.mathworks.com/help/comm/ug/error-detection-and-correction.html
%* 
%* Author: T.J. Cheng, 2016
%* 
%*   2016-11-20: Error: ldpc_enc_IRA(s, N, rate, q1)
%*   2016-11-21: Bug fixed. Everything works fine so far.
%*   2016-11-25: Remove built-in encoder (not compatible to QC-LDPC)
%*   2016-11-26: Adopts the method of syndrome-check described in [2]
%*

clc
clear all;

% params
N       = 64800;
rate    = 9/15;
q1      = 72;

% generate information bits
s = randi([0 1], 1, N * rate);
%disp(s);

% generate parity-check matrix
disp('generating parity-check matrix...');
H = ldpc_pcmg_B(N, rate, q1);
fprintf('\tdone\n');
spy(H);

% encode information bits
c = ldpc_enc_B(s, N, rate, q1);

% perform syndrome-check
disp('performing syndrome-check...');
fprintf('\tdone with result = %d\n', sum(rem(c * H', 2)));