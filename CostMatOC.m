function B = CostMatOC(A, isON)
% Function: construct weighed bipartite network 
% INPUT:
% A: adjacent matrix: %aij: i-j
% isON: is output nodes???
% OUTPUT:
% B: weighted bipartite graph

% Reference
% Wu, L., Shen, Y., Li, M., Wu, F.X.: Network output controllability-based method
% for drug target identification. IEEE Trans Nanobioscience 14(2), 184¨C191 (2015)

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
sigma = 1+ 0.1/num;
epsilon = 0.1/num;

for i = 1:num
    for j = 1:num
        aij = A(i, j); %j->i
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

wru = zeros(num, 1);
wru(isON==1) = sigma;
wru(isON==0) = epsilon;

Wru = diag(wru);
B = [B, Wru];

end