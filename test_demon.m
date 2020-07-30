clear
clc
close

%% 
% A(matrix): where aij:i->j
% isON(vec): output nodes
% isPN(vec): preference nodes
% isCN(vec): constrained nodes
% Drugs(matrix): each column represent one 
% DrugsW(vec): weights of drugs
load test_network1

% MinISOC = CtrMinISOC(A, isON, drugs)
COC = CtrCOC(A, isON, isCN)






