function B = CostMatPCOC(A, isON, CN)
% FUNCTION: Generate weighted bipartite network for PCOC
% INPUT:
% A: adjacent matrix: %aij: i->j
% isON: is output nodes
% CN: weighted constrained nodes
% OUTPUT:
% B: weighted bipartite network

% Reference
% Hao Gao*, Xiangmao Meng*, Min Li and Fang-Xiang Wu, MinISOC: Minimum Input 
% Selections Output Control Problems for Drug Combinations Identification
% (In preparation)

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
epsilon = 1+CN;

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

idx = find(CN);
Wru = zeros(num, length(idx));
for i = 1:length(idx)
    Wru(idx(i), i) = epsilon(i);
end

B = [B, Wru];

end