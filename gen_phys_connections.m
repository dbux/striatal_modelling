function[connections, list] = gen_phys_connections(striatum, attr)
    %  Humphries, Wood & Gurney (2010)
    %  Page 10

    warning('off', 'stats:pdist2:DataConversion');
        
    % Extract necessary attributes
    phys = attr.phys;
    conn = attr.conn;
    flags = attr.flags;

    reverseStr = '';            % Required string for progress indicator
    E = gen_e;                  % Generate values for E
    ptr = 1;                    % Points to the next available synaptic connection entry   
    
    msn = 1;                    % How neurons will be represented numerically            
    fsi = 3;

    % Prellocate synaptic connection lists
    connections.msnmsn = zeros((1000 * length(striatum.linear)), 3);
    connections.fsimsn = zeros((1000 * length(striatum.linear)), 3);
    connections.fsifsi = zeros((1000 * length(striatum.linear)), 3);
    connections.gap    = [];    % Gap connection list does not need preallocation

    if flags.progress
        timer.conn = tic;
        fprintf('\nCreating connections…')
        fprintf('\nFor a full-size network (1mm³) this may take some time')
        fprintf('\nStart time: %s | ', datestr(now, 'HH:MM:SS'))
    end

    % For each neuron…
    for i = 1:length(striatum.linear)
        timer.neuron = tic;

        % Connection information
        coninf.count    = 0;    % Counts number of connections created 
        coninf.dist_all = [];	% Distance from this neuron to all other neurons
        coninf.exp_all  = [];   % Expected number of connections of all types to all other neurons 
        coninf.exp_this = [];   % Expected number of connections of appropriate types to all other neurons
        coninf.trials   = [];   % Number of connection trials to perform for all potential connections
        coninf.prob     = [];   % Probability of making a connection at each trial to all other neurons

        if flags.debug
            fprintf('\nConnecting neuron %d of %d… ', i, length(striatum.linear))
        end

        % Get the distance from current neuron to all other neurons. Must be a
        % round number to look up expected connections from E later.    
    %     coninf.dist_all = round(distancePoints3d(striatum.neurons(i,:), striatum.neurons));
        coninf.dist_all = pdist2(striatum.neurons(i,:), striatum.neurons);       

        % This is a bit of a hack - to prevent errors from unwanted 0-distance
        % connections (i.e. self-connections) we set the distance to 2000
        % (maximum distance, with 0 expected connections for all types)
        coninf.dist_all(coninf.dist_all==0) = length(E);

        % Get the type of connection to all other neurons
        % MSN = 1, FSI = 3
        % So MSN-MSN = 1+1 = 2; FSI-MSN = 3+1 = 4; FSI-FSI = 3+3 = 6
        % (Must be converted to UINT32 to avoid forcing connection lists into UINT8 format)
        coninf.type = uint32(striatum.linear(i) + striatum.linear);

        % If current neuron is 1 (MSN), remove connections of type 4 (MSN-FSI) as they are invalid
        if striatum.linear(i) == msn
            coninf.type(coninf.type==(msn+fsi)) = 0;
        end 

        % Get the expected number of connections from this neuron to every
        % other neuron based on their relative distance. At this stage we don't
        % discriminate based on the type of connection being made, and obtain
        % expected connection values for all connection types.
        coninf.exp_all = E(:,round(coninf.dist_all))';

        % Get the values of only the desired connection types. For FSI-FSI 
        % connections we want values for both synapses and gap junctions so two 
        % columns are needed
        coninf.exp_this(coninf.type==(msn+msn), 1) = coninf.exp_all(coninf.type==(msn+msn), 1);   
        coninf.exp_this(coninf.type==(fsi+msn), 1) = coninf.exp_all(coninf.type==(fsi+msn), 2);   
        coninf.exp_this(coninf.type==(fsi+fsi), 1) = coninf.exp_all(coninf.type==(fsi+fsi), 3);   
        coninf.exp_this(coninf.type==(fsi+fsi), 2) = coninf.exp_all(coninf.type==(fsi+fsi), 4);

        % Because some expected connection values are >1, a standard method for
        % generating probabilistic connections is required:
        % Let n be the expected number of connections between two neurons.
        % (2*ceiling-n) trials are run, each with a (0.5 * n/ceiling-n) chance 
        % of making a connection.
        coninf.trials = 2 * ceil(coninf.exp_this);
        coninf.prob = 0.5 * (coninf.exp_this ./ ceil(coninf.exp_this));

        % Run the trials for each pair and generate connections. Store the
        % connection pair, type, and the distance between the two neurons
        
        % for (length of [trials])
        % generate trials(n) number of rand numbers
        % if any of them are lower than prob(n)
        % create that many connections between i and (target neuron)
        % Apply all that to every element in trials with arrayfun?
%         for n = 1: numel(coninf.trials)
%         if any(rand(1, coninf.trials) < coninf.prob

        for k = 1:length(coninf.trials)
            for j = 1:coninf.trials(k)
                if rand < coninf.prob(k, 1)
                    coninf.dist_this = coninf.dist_all(k);
                    if coninf.type(k) == (msn+msn)
                        connections.msnmsn(ptr,:) = [i, k, coninf.dist_this];
                    elseif coninf.type(k) == (fsi+msn)
                        connections.fsimsn(ptr,:) = [i, k, coninf.dist_this];
                    elseif coninf.type(k) == (fsi+fsi)
                        connections.fsifsi(ptr,:) = [i, k, coninf.dist_this];
                    else
                        error('Unexpected connection type!')
                    end
                    ptr = ptr+1;
                    coninf.count = coninf.count + 1;
                end
            end
        end

        % If the current neuron is an FSI, run connection trials for gap junctions
        if striatum.linear(i) == fsi
            for j = 1:max(coninf.trials(:,2))
                for k = 1:length(coninf.trials)
                    if coninf.trials(k,2) ~= 0
                        if rand < coninf.prob(k,2)
                            coninf.dist_this = coninf.dist_all(k);
                            connections.gap = [connections.gap; i k coninf.dist_this];
                            coninf.count = coninf.count + 1;
                        end
                        coninf.trials(k,2) = coninf.trials(k,2) - 1;
                    end
                end
            end
        end

        if flags.debug
            fprintf('made %d connections in %1.2f seconds', coninf.count, toc(timer.neuron))
        end

        if flags.debug == 0 && flags.progress
            percentDone = 100 * i / length(striatum.linear);
            msg = sprintf('Percent done: %3.3f', percentDone);
            fprintf([reverseStr, msg])
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end

    % Remove empty space at the end of connection lists
    connections.msnmsn = connections.msnmsn(any(connections.msnmsn,2),:);
    connections.fsimsn = connections.fsimsn(any(connections.fsimsn,2),:);
    connections.fsifsi = connections.fsifsi(any(connections.fsifsi,2),:);

    % Convert connections distances to delay
    % Multiply conduction velocity in m/s by 1,000 for equivalent in μm/ms
% 	connections.msnmsn(:,3) = max(connections.msnmsn(:,3) ./ (phys.cv_msnmsn * 1000), conn.delay_min);
% 	connections.fsimsn(:,3) = max(connections.fsimsn(:,3) ./ (phys.cv_fsimsn * 1000), conn.delay_min);
% 	connections.fsifsi(:,3) = max(connections.fsifsi(:,3) ./ (phys.cv_fsifsi * 1000), conn.delay_min);
    
    % Extract list of each type of neuron
    [connections, list] = get_neuron_list(connections); 

    if flags.progress
        fprintf('\nCreating connections took %1.2f minutes\n', toc(timer.conn)/60 )
    end

%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     %%%%%% END FUNCTION HERE AND MOVE BELOW ELSEWHERE
% 
%     % % % % % % % % %
%     % Generate connectivity statistics
%     statname = [striatum.dirname '/connection_stats.csv'];
%     fc = fopen(statname, 'at+');
%     fprintf(fc, '\n%d MSNs, %d FSIs, FSI ratio %d%%', ...
%         length(find(striatum.linear==msn)), length(find(striatum.linear==fsi)), phys.fsi_pct);
%     fprintf(fc, '\nConnection type, No. contacts, STD, Distance (µm), STD');
% 
%     % MSNs - 1 MSN
%     try
%         msn1msn = gen_conn_stats(connections.msnmsn, striatum, 0);   
%         fprintf(fc, '\nMSNs - 1 MSN, %1.2f, %1.2f, %1.2f, %1.2f', mean(msn1msn.numbers), std(msn1msn.numbers), ...
%             mean(msn1msn.dists), std(msn1msn.dists));
%     catch
%         fprintf(fc, '\nError generating MSNs - 1 MSN connection statistics!');
%     end
% 
%     % FSIs - 1 MSN
%     try
%         fsi1msn = gen_conn_stats(connections.fsimsn, striatum, 0);
%         fprintf(fc, '\nFSIs - 1 MSN, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsi1msn.numbers), std(fsi1msn.numbers), ...
%             mean(fsi1msn.dists), std(fsi1msn.dists));
%     catch
%         fprintf(fc, '\nError generating FSIs - 1 MSN connection statistics!');
%     end
% 
%     % 1 FSI - MSNs
%     try
%         fsimsns = gen_conn_stats(connections.fsimsn, striatum, 1);
%         fprintf(fc, '\n1 FSI - MSNs, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsimsns.numbers), std(fsimsns.numbers), ...
%             mean(fsimsns.dists), std(fsimsns.dists));
%     catch
%         fprintf(fc, '\nError generating 1 FSI - MSNs connection statistics!');
%     end
% 
%     % FSIs - 1 FSI
%     try
%         fsi1fsi = gen_conn_stats(connections.fsifsi, striatum, 0);
%         fprintf(fc, '\nFSIs - 1 FSI, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsi1fsi.numbers), std(fsi1fsi.numbers), ...
%             mean(fsi1fsi.dists), std(fsi1fsi.dists));
%     catch
%         fprintf(fc, '\nError generating FSIs - 1 FSI connection statistics!');
%     end
% 
%     % FSI gap
%     try
%         fsigap = gen_conn_stats(connections.gap, striatum, 0);
%         fprintf(fc, '\nFSI gap, %1.2f, %1.2f, %1.2f, %1.2f\n', mean(fsigap.numbers), std(fsigap.numbers), ...
%             mean(fsigap.dists), std(fsigap.dists));
%     catch
%         fprintf(fc, '\nError generating FSI gap connection statistics!');
%     end
% 
%     fclose(fc);
%    
% 
%     % DISTANCES ARE WRONG
% 
%     if flags.progress
%         fprintf('\nStriatal connectivity statistics (%d%% FSI):', phys.fsi_pct)
%         fprintf('\n----------------------------------------------------------')
%         fprintf('\nType\t\t| No. of contacts\t| Distance (µm)')
%         fprintf('\n----------------------------------------------------------\n')
%         fprintf('MSNs - 1 MSN\t| %1.2f, ± %1.2f\t| %1.2f ± %1.2f\n', mean(msn1msn.numbers), std(msn1msn.numbers), ...
%             mean(msn1msn.dists), std(msn1msn.dists))
%         fprintf('FSIs - 1 MSN\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsi1msn.numbers), std(fsi1msn.numbers), ...
%             mean(fsi1msn.dists), std(fsi1msn.dists))
%         fprintf('1 FSI - MSNs\t| %1.2f, ± %1.2f\t| %1.2f ± %1.2f\n', mean(fsimsns.numbers), std(fsimsns.numbers), ...
%             mean(fsimsns.dists), std(fsimsns.dists))
%         fprintf('FSIs - 1 FSI\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsi1fsi.numbers), std(fsi1fsi.numbers), ...
%             mean(fsi1fsi.dists), std(fsi1fsi.dists))
%         try
%             fprintf('FSI gap junc\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsigap.numbers), std(fsigap.numbers), ...
%             mean(fsigap.dists), std(fsigap.dists))
%         catch
%         end
%         fprintf('----------------------------------------------------------\n')
%         fprintf('Data saved to file. Have a nice day!\n')
%     end
%     

    %% FUNCTIONS
    function E = gen_e
        % Generates the value E for each connection type as specified on page 9 of
        % Humphries, Wood & Gurney (2010)

        % Set maximum distance to greater than striatal edge length to allow for 
        % corner-to-corner checks
        d = 1:2000;
        
        % Values from page 10, table 4
        alpha = [0.511, -0.921, -0.695, 1.322];
        beta  = [1.033, 1.033, 1.38, 2.4];
        gamma = [0.042, 0.042, 0.057, 0.016];
        delta = [26.8, 26.8, 15.6, 43.3];
        eta   = [0.0039, 0.0039, 0.0036, 0.0029];
        
        % Preallocate
        E = zeros(length(alpha), length(d));

        % Page 9, equation 20
        for c = 1:numel(alpha) 
            E(c,:) = exp(-alpha(c) - beta(c) .* (1 - exp(-gamma(c) .* (d - delta(c)))) .* exp(eta(c) .* d));
        end
    end
end