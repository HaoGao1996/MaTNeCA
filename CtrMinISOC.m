function MinISOC = CtrMinISOC(A, isON, drugs, w)
% FUNCTION: Minimum Input Selection Output Control
% INPUT:
% A: the directed network, where aij: i->j
% isON: is output nodes?
% drugs: each column represents a drug, each row indicates a node targeted
% by drugs
% w: weights of drugs
% OUTPUT:
% MinISPOC: the set of drugs

% Reference
% Hao Gao*, Xiangmao Meng*, Min Li and Fang-Xiang Wu, MinISOC: Minimum Input 
% Selections Output Control Problems for Drugs Indentification
% (In preparation)

% Hao Gao. MaTNeCA: MatLab Toolkit for Network Control Analysis. 2020

% Copyright: Hao Gao (Hougogh)
% Contact: ggogh111@gmail.com
% Date: 2020/08/03

% ---version 0.11---

%%
A = full(A); %aij:i->j
num_n = length(A);
num_d = size(drugs, 2); 
A = [A, zeros(num_n, num_d); drugs', zeros(num_d)];
isON = [isON; zeros(num_d, 1)];

if nargin < 3
    print('Error!!!! Wrong Number of Arguments!!')
    return
elseif nargin == 3
    % if no w as inputs, they have the same weights.
    w = ones(num_d, 1);
end

w = 0.8 * w / (sum(w) * num_n);% normalization

% the drugs (nodes) are considered as the constrained set of nodes
CN = [zeros(num_n, 1); w]; 

% get the preferential constratrained output control results
PCOC = CtrPCOC(A, isON, CN);

% get the label of drugs
MinISOC = PCOC - num_n;

end