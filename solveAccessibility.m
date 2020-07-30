function sf = solveAccessibility(A, so, isON, isCN)
% FUNCTION: solve accessibility of output control
% INPUT:
% A: is the directed network, where aij: i->j
% so: original steering nodes
% isON: is output nodes
% isCN: is constrained nodes
% OUTPUT:
% sf: final steering nodes

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/07/30

%% if the input nodes have no constrains
if narin == 3
    A = sparse(A);
    inAccessibleNode = get_inAccessibleNode(A, so, isON);
    
    if ~isempty(inAccessibleNode)
        c = components(A);
        c_ON = []; % if one node is accessible, all nodes are accessible
        st = []; % temporal steering nodes
        for i = 1:length(inAccessibleNode)
            if isempty(find(c_ON == c(inAccessibleNode(i))))
                c_ON = [c_ON; c(inAccessibleNode(i))];
                
                % if no constraint, choose any node
                st = [st; inAccessibleNode(i)]; 
            end
        end
        sf = [so; st];
    else
        sf = so;
    end
%% if the input nodes have constrains   
else
    A = sparse(A);
    inAccessibleNode = get_inAccessibleNode(A, so, isON);
    sf = so;
    
    while(~isempty(inAccessibleNode))
        % Identify the constained nodes that can reach inAccessibleNode(1)
        d = shortest_paths(A', inAccessibleNode(1));
        d(d~=inf) = 0; % inaccessible
        d(d==inf) = 1; % accessible
        isAccessibility = d.*isCN;
        inAccessibleCN = find(isAccessibility);
        sf = [sf; inAccessibleCN(1)];
        
        inAccessibleNode = get_inAccessibleNode(A, sf, isON);     
        
    end
    
end

end
