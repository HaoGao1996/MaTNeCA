function COC = CtrCOC(A, isON, isCN)
% FUNCTION: constrained output control
% INPUT:
% A: is the directed network, where aij: i->j
% isON: is output nodes
% isCN: is constrained nodes
% OUTPUT:
% COC

% Reference:
% Guo, W.F., Zhang, S.W., Wei, Z.G., Zeng, T., Liu, F., Zhang, J., Wu, F.X., Chen,
% L.: Constrained target controllability of complex networks. J. Stat. Mech. Theory
% Exp. 2017(6), 063402 (2017)

% Hao Gao. MaTNeCA: MatLab Toolkit for Network Control Analysis. 2020

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/07/30

%% Maximum Weighted Complete Matching (MWCM)
A = full(A);
B = CostMatCOC(A, isON, isCN);


[col4row,row4col,gain,u,v] = assign2D(B, true);

m = row4col(length(A)+1:end);
m(m==0) = [];
COC1 = m;

%% Accessibility test from constrained nodes
COC = solveAccessibility(A, COC1, isON, isCN);

end