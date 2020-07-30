function inAccessibleNode = get_inAccessibleNode(A, so, isON)
% FUNCTION: get inaccessible nodes for a given set of input nodes and 
%           target nodes
% INPUT:
% A: is the directed network, where aij: i->j
% so: original steering nodes
% isON: is output nodes
% OUTPUT:
% inAccessibleNode: inaccessible nodes

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/07/30

%%
d = [];
for i=1:length(so)
    di = shortest_paths(A, so(i));
    d = [d, di];
end

[d, ~] = min(d, [], 2);

d(d~=inf) = 0; % inaccessible
d(d==inf) = 1; % accessible

isAccessibility = d.*isON;
inAccessibleNode = find(isAccessibility);

end