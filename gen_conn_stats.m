function[output] = gen_conn_stats(conn, striatum, bool)
% The following code generates statistics similar to 
% Humphries, Wood & Gurney (2010), page 11, table 5. It operates as follows:

% x.dists = Get the distance (col 3) from the appropriate connection list
% for every neuron that is in the centre region. Convert the result from
% UINT to single to avoid problems.

% x.list = Get a list of all the times a neuron exists as either an
% afferent (col 1) or efferent (col 2) for all neurons in the centre region.
% E.g for MSNs - 1 MSN, will list every MSN in the centre region innervated
% by at least 1 other MSN (located anywhere), repeated by as many incoming
% MSN connections exist.

% x.members = As x.list but with duplicates removed
% x.numbers = Find out how many contacts exist for each neuron in x.members

% Boolean value indicates if starting or destination population is single

if bool
    % For 1 x - multiple y
    output.dists = single(conn(find(ismember(conn(:,1), striatum.linear_centre(:,1)))',3));
    output.list = single(conn(find(ismember(conn(:,1), striatum.linear_centre(:,1)))',1));
else
    % For multiple x - 1 y
    output.dists = single(conn(find(ismember(conn(:,2), striatum.linear_centre(:,1)))',3));
    output.list = single(conn(find(ismember(conn(:,2), striatum.linear_centre(:,1)))',2));
end

output.members = unique(output.list);
output.numbers = hist(output.list, output.members);
output.numbers = output.numbers(output.numbers ~= 0);