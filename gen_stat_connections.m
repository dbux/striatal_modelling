function[connections, list] = gen_stat_connections(attr)
% Generates statistically-based striatal connections

    % Extract necessary attributes
    stat = attr.stat;
    conn = attr.conn;
    progress = attr.flags.progress;

    if progress
        timer.stat = tic;
        fprintf('\nCreating statistical connections…')      
    end

    % Given the number of neurons and the number of expected connections, what
    % is the probability of a connection between any two neurons?
    prob.msnmsn = stat.con_msnmsn / stat.num_msn;
    prob.fsimsn = stat.con_fsimsn / stat.num_msn;
    prob.fsifsi = stat.con_fsifsi / stat.num_fsi;
    prob.fsigap = stat.con_fsigap / stat.num_fsi;
    
    % Generate connection (note gap junctions are a subset of FSI GABA connections)
    connections.msnmsn = generate_connections(stat.num_msn, stat.num_msn, prob.msnmsn);
    connections.fsimsn = generate_connections(stat.num_fsi, stat.num_msn, prob.fsimsn);
    connections.fsifsi = generate_connections(stat.num_fsi, stat.num_fsi, prob.fsifsi);
    connections.gap    = generate_connections([], [], prob.fsigap, connections.fsifsi);
    
    % Extract list of each type of neuron
    [connections, list] = get_neuron_list(connections); 
    
    if progress
        fprintf('\nCreating connections took %1.2f seconds\n', toc(timer.stat))
    end

    %% FUNCTIONS
    function connections = generate_connections(num_source, num_dest, prob, varargin) 
        % Create connections based on number of source and destination neurons
        % and their connections probability
        
        if numel(varargin) > 0 
            % If a connection list is supplied, use this as the starting point
            connections = varargin{1};
        else
            % Create an all-to-all connection list between source and destination neurons
            connections = gen_list_allall(num_source, num_dest);
        end
            
        % Append random connection chance
        connections = add_random_chance(connections);

        % If connection probability is greater than one (i.e. more than one
        % connection should be made) then connections with a probability value
        % less than (probability - 1) are duplicated and none are removed
        if prob > 1
            prob_curr = prob;
            while prob_curr > 1
                conn_extra = connections((connections(:, 3) <= prob_curr - 1), :);
                connections = [connections ; conn_extra];
                prob_curr = prob_curr - 1;
            end
        else
            % Trim the connections to match the expected probability value
            connections = connections((connections(:, 3) <= prob), :);
        end
        
        % Remove connection probability
        connections(:, 3) = [];

%         if numel(varargin) > 0
%             % If a second argument is supplied use this as the connection delay,
%             % otherwise set to null (as required for gap connections)
%             try
%                 connections(:, 3) = varargin{2};
%             catch
%                 connections(:, 3) = [];
%             end
%         else
%             % Replace the random probability value with a fixed connection delay
%             connections(:, 3) = conn.delay_min;
%         end
        
        
        function[list] = gen_list_allall(n1, n2)
            % Create all-to-all connection lists
            if length(n1) == 1
                n1 = 1:n1;
            end

            if length(n2) == 1
                n2 = 1:n2;
            end

            % Get length of each array
            len1 = length(n1);
            len2 = length(n2);

            % Create repeating lists containing all values the required number of times
            l1 = repmat(n1', len2, 1);
            l2 = repmat(n2', 1, len1)';
            l2 = l2(:);

            % Create all-to all list
            list = [l1' ; l2']';
        end
        
 
        function list = add_random_chance(list)
            % Create random connection chance for every neuron pair
            % (Must be done in two steps or every connection gets the same value)
            r = rand(length(list), 1); 
            list(:,3) = r;
        end
    end
end


  