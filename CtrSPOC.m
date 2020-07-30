function SPOC = CtrSPOC(A, isON, isPN)
% FUNCTION: set preference output control
% INPUT:
% A: is the directed network, where aij: i->j
% isON: is output nodes
% isPN: is preferential nodes
% OUTPUT:
% SPOC: 

% Reference:
% Hao Gao, Min Li and Fang-Xiang Wu. SPOC: Identification of Drug Targets in 
% Biological Networks via Set Preference Output Control. ISBRA2020 (Accepted)

% Hao Gao. MaTNeCA: MatLab Toolkit for Network Control Analysis. 2020

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/07/30

%% Maximum Weighted Complete Matching (MWCM)
A = full(A);
B = CostMatSPOC(A, isON, isPN);


[col4row,row4col,gain,u,v] = assign2D(B, true);

m = row4col(length(row4col)/2+1:end);
m(m==0) = [];
SPOC1 = m;

%% accessibility
SPOC = solveAccessibility(A, SPOC1, isON);

end