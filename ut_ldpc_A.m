%* Unit Test for Type A (MET-Structured) LDPC Encoder
%*
%* References:
%*   [1] Physical Layer Protocol, document ATSC A/322, ATSC, Sep. 2016
%*   [2] MathWorks, "Error Detection and Correction", Available: 
%*         https://www.mathworks.com/help/comm/ug/error-detection-and-correction.html
%* 
%* Author: T.J. Cheng, 2016
%* 
%*   2016-11-26: Everything works fine so far.
%*   2016-11-26: Adopts the method of syndrome-check described in [2]
%*

clc
clear all;

% params
N       = 64800;
rate    = 7/15;
M1      = 1080;
M2      = 33480;
q1      = 3;
q2      = 93;

% generate information bits
s = randi([0 1], 1, N * rate);
%disp(s);

% generate parity-check matrix
disp('generating parity-check matrix...');
H = ldpc_pcmg_A(N, rate, M1, M2, q1, q2);
% H2 = ldpc_pcmg_MET2(N, rate, M1, M2, q1, q2);
fprintf('\tdone\n');
% isequal(H, H2)
spy(H);

% encode information bits
% c = ldpc_enc_MET(s, N, rate, M1, M2, q1, q2);
c = ldpc_enc_A(s, N, rate, M1, M2, q1, q2);

% perform syndrome-check
disp('performing syndrome-check...');
fprintf('\tdone with result = %d\n', sum(rem(c * H', 2)));