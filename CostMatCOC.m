function B = CostMatCOC(A, isON, isCN)
% FUNCTION: Generate weighted bipartite network
% INPUT:
% A: adjacent matrix: %aij: i->j
% isON: is output nodes???
% isCN: is constained nodes???
% OUTPUT:
% B: weighted bipartite network

% Reference:
% Guo, W.F., Zhang, S.W., Wei, Z.G., Zeng, T., Liu, F., Zhang, J., Wu, F.X., Chen,
% L.: Constrained target controllability of complex networks. J. Stat. Mech. Theory
% Exp. 2017(6), 063402 (2017)

% Hao Gao. MaTNeCA: MatLab Toolkit for Network Control Analysis. 2020

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.coms
% Date: 2020/08/03

% ---version 0.11---

%%
A = A'; % aij j->i
num = length(A);
B = zeros(num, num);

alpha = 2;
beta = 1;
epsilon = 1+0.1/num;

for i = 1:num
    for j = 1:num
        aij = A(i, j);
        if isON(i)
           if aij
               B(i, j) = alpha;
           end
        else
           if aij
               B(i, j) = beta;
           else
               if i==j
                   B(i, j) = beta;
               end
           end
        end
    end
end

idx = find(isCN);
Wru = zeros(num, length(idx));
for i = 1:length(idx)
    Wru(idx(i), i) = epsilon;
end

B = [B, Wru];

end