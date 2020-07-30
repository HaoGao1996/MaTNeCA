function PCOC = CtrPCOC(A, isON, CN)
% FUNCTION: constrained output control
% INPUT:
% A: is the directed network, where aij: i->j
% isON: is output nodes
% CN: weighted constrained nodes
% OUTPUT:
% PCOC

% Reference
% Hao Gao*, Xiangmao Meng*, Min Li and Fang-Xiang Wu, MinISOC: Minimum Input 
% Selections Output Control Problems for Drug Combinations Identification
% (In preparation)

% Hao Gao. MaTNeCA: MatLab Toolkit for Network Control Analysis. 2020

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/07/30

%% Maximum Weighted Complete Matching (MWCM)
A = full(A);
B = CostMatPCOC(A, isON, CN);


[col4row,row4col,gain,u,v] = assign2D(B, true);

m = row4col(length(A)+1:end);
m(m==0) = [];
PCOC1 = m;

%% Accessibility
PCOC = solveAccessibility(A, PCOC1, isON, CN);

end