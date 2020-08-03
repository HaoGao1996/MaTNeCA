function B = CostMatSPOC(A, isON, isPN)
% FUNCTION: Generate weighted bipartite network
% INPUT:
% A: adjacent matrix: %aij: i->j
% isON: is output nodes???
% isPN: is preferential nodes???
% OUTPUT:
% B: weighted bipartite network

% Reference:
% Hao Gao, Min Li and Fang-Xiang Wu. SPOC: Identification of Drug Targets in 
% Biological Networks via Set Preference Output Control. ISBRA2020 (Accepted)

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
epsilon = 1+0.5/num;
gamma = 1+0.1/num;


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

wru = zeros(num, 1);
for i = 1:num
	if isPN(i) && isON(i)
		wru(i) = epsilon;
	end

	if isPN(i) && ~isON(i)
		wru(i) = epsilon-1;
	end

	if ~isPN(i) && isON(i)
		wru(i) = gamma;
	end
	
	if ~isPN(i) && ~isON(i)
		wru(i) = gamma-1;
	end
end

Wru = diag(wru);
B = [B, Wru];

end