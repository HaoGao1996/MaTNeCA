clear
clc
close
%%
% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/08/03

% ---version 0.11---

%% 
% A(matrix): where aij:i->j
% isON(vec): output nodes
% isPN(vec): preference nodes
% isCN(vec): constrained nodes
% drugs(matrix):  each column represents a drug, each row indicates a node targeted
% by drugs 
% w(vec): weights of drugs

%%
load test_network1

MinISOC = CtrMinISOC(A, isON, drugs)
% COC = CtrCOC(A, isON, isCN)






