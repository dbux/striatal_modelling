function[list] = gen_list_allall(n1, n2)
% If inputs are integers instead of lists, create lists from 0 to (input
% value - 1) to account for SpineCreator indexing
if length(n1) == 1
    n1 = 0:n1 - 1;
end

if length(n2) == 1
    n2 = 0:n2 - 1;
end

% Get length of each array
len1 = length(n1);
len2 = length(n2);

% Create repeating lists containing all values the required number of times
l1 = repmat(n1', len2, 1);
l2 = repmat(n2', 1, len1)';
l2= l2(:);

% Create all-to all list
list = [l1' ; l2']';