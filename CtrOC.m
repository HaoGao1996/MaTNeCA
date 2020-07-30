function OC = CtrOC(A, isON)
% FUNCTION: Output control
% According to the article published by Wu et al.
% INPUT:
% A: is the directed network, where aij: i->j
% isON: is output nodes
% Output:
% OC: the driver nodes of output control

% Reference
% Wu, L., Shen, Y., Li, M., Wu, F.X.: Network output controllability-based method
% for drug target identification. IEEE Trans Nanobioscience 14(2), 184¨C191 (2015)

% Hao Gao. MaTNeCA: MatLab Toolkit for Network Control Analysis. 2020

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/07/30

%% Maximum Weighted Complete Matching (MWCM)
% dilation
A = full(A);
B = CostMatOC(A, isON);

[col4row,row4col,gain,u,v]=assign2D(B, true);

m = row4col(length(row4col)/2+1:end);
m(m==0) = [];
OC1 = m;

%% Accessibility

OC = solveAccessibility(A, OC1, isON);

end