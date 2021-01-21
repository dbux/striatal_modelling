function gen_phys_connstats(striatum, connections)

    % MSNs - 1 MSN
    try
        msn1msn = get_stats(connections.msnmsn, striatum, 0);   
%         fprintf(fc, '\nMSNs - 1 MSN, %1.2f, %1.2f, %1.2f, %1.2f', mean(msn1msn.numbers), std(msn1msn.numbers), ...
%             mean(msn1msn.dists), std(msn1msn.dists));
    catch
        fprintf('\nError generating MSNs - 1 MSN connection statistics!');
    end

%     connstats = table;

%     % % % % % % % % %
%     % Generate connectivity statistics
%     statname = [striatum.dirname '/connection_stats.csv'];
%     fc = fopen(statname, 'at+');
%     fprintf(fc, '\n%d MSNs, %d FSIs, FSI ratio %d%%', ...
%         length(find(striatum.linear==msn)), length(find(striatum.linear==fsi)), phys.fsi_pct);
%     fprintf(fc, '\nConnection type, No. contacts, STD, Distance (µm), STD');



    % FSIs - 1 MSN
    try
        fsi1msn = get_stats(connections.fsimsn, striatum, 0);
        fprintf(fc, '\nFSIs - 1 MSN, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsi1msn.numbers), std(fsi1msn.numbers), ...
            mean(fsi1msn.dists), std(fsi1msn.dists));
    catch
        fprintf(fc, '\nError generating FSIs - 1 MSN connection statistics!');
    end

    % 1 FSI - MSNs
    try
        fsimsns = gen_conn_stats(connections.fsimsn, striatum, 1);
        fprintf(fc, '\n1 FSI - MSNs, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsimsns.numbers), std(fsimsns.numbers), ...
            mean(fsimsns.dists), std(fsimsns.dists));
    catch
        fprintf(fc, '\nError generating 1 FSI - MSNs connection statistics!');
    end

    % FSIs - 1 FSI
    try
        fsi1fsi = gen_conn_stats(connections.fsifsi, striatum, 0);
        fprintf(fc, '\nFSIs - 1 FSI, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsi1fsi.numbers), std(fsi1fsi.numbers), ...
            mean(fsi1fsi.dists), std(fsi1fsi.dists));
    catch
        fprintf(fc, '\nError generating FSIs - 1 FSI connection statistics!');
    end

    % FSI gap
    try
        fsigap = gen_conn_stats(connections.gap, striatum, 0);
        fprintf(fc, '\nFSI gap, %1.2f, %1.2f, %1.2f, %1.2f\n', mean(fsigap.numbers), std(fsigap.numbers), ...
            mean(fsigap.dists), std(fsigap.dists));
    catch
        fprintf(fc, '\nError generating FSI gap connection statistics!');
    end

    fclose(fc);
   

    % DISTANCES ARE WRONG

    if flags.progress
        fprintf('\nStriatal connectivity statistics (%d%% FSI):', phys.fsi_pct)
        fprintf('\n----------------------------------------------------------')
        fprintf('\nType\t\t| No. of contacts\t| Distance (µm)')
        fprintf('\n----------------------------------------------------------\n')
        fprintf('MSNs - 1 MSN\t| %1.2f, ± %1.2f\t| %1.2f ± %1.2f\n', mean(msn1msn.numbers), std(msn1msn.numbers), ...
            mean(msn1msn.dists), std(msn1msn.dists))
        fprintf('FSIs - 1 MSN\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsi1msn.numbers), std(fsi1msn.numbers), ...
            mean(fsi1msn.dists), std(fsi1msn.dists))
        fprintf('1 FSI - MSNs\t| %1.2f, ± %1.2f\t| %1.2f ± %1.2f\n', mean(fsimsns.numbers), std(fsimsns.numbers), ...
            mean(fsimsns.dists), std(fsimsns.dists))
        fprintf('FSIs - 1 FSI\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsi1fsi.numbers), std(fsi1fsi.numbers), ...
            mean(fsi1fsi.dists), std(fsi1fsi.dists))
        try
            fprintf('FSI gap junc\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsigap.numbers), std(fsigap.numbers), ...
            mean(fsigap.dists), std(fsigap.dists))
        catch
        end
        fprintf('----------------------------------------------------------\n')
        fprintf('Data saved to file. Have a nice day!\n')
    end

    function[output] = get_stats(conn, striatum, bool)
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
            output.dists = single(conn(find(ismember(conn(:, 1), striatum.linear_centre(:, 1)))', 3));
            output.list  = single(conn(find(ismember(conn(:, 1), striatum.linear_centre(:, 1)))', 1));
        else
            % For multiple x - 1 y
            output.dists = single(conn(find(ismember(conn(:, 2), striatum.linear_centre(:, 1)))', 3));
            output.list  = single(conn(find(ismember(conn(:, 2), striatum.linear_centre(:, 1)))', 2));
        end

        output.members = unique(output.list);
        output.numbers = hist(output.list, output.members);
        output.numbers = output.numbers(output.numbers ~= 0);
    end
end