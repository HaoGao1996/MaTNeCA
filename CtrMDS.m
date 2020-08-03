function MDS = CtrMDS(A)
% FUNCTION: return the matching of associated bipartite networks
% Input: 
% A is the directed network, where aij: i->j; !!!!!!!!!!!!!
% Output: 
% MDS: minimum set of driver nodes

% Reference:
% Liu, Y.Y., Slotine, J.J., Barab?asi, A.L.: Controllability of complex 
% networks. Nature 473(7346), 167 (2011)

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.coms
% Date: 2020/08/03

% ---version 0.11---

%% Construct bipartite graph
A = A'; % Aij: j->i

num = length(A);
B = zeros(2 * num);
A = full(A);

for i = 1:num
    for j = 1:num
        if A(i, j)
            % symmetric
            B(i, j+num) = 1;
            B(j+num, i) = 1;
        end
    end
end
B = sparse(B);

%% Maximum Matching
m0 = edmonds_maximum_cardinality_matching(B);
m = m0(1:num);
m(find(m==0)) = num;
m = m - num;

% Unmatched nodes are driver nodes
MDS = find(m==0); 

end