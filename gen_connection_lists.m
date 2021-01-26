function[output_lists, attr] = gen_connection_lists(connections, list, attr, varargin)
   
    % Extract necessary attributes
    phys = attr.phys;
    conn = attr.conn;
    flags = attr.flags;
    
    % Path for saving lists (and checking existence of current lists)
    attr.list_path = fullfile(attr.save_path, 'connection_lists', sprintf('delay_%1.1f', conn.delay_mult));

	if attr.flags.physical
        try
            striatum = varargin{1};
        catch
            error('No striatum found - cannot generate connection lists') 
        end
	end
    
    % Generate connection ID prefix
    conn_prefix = matlab.lang.makeValidName(strcat('bkMSN', num2str(conn.bkg_msn), '_bkFSI', num2str(conn.bkg_fsi)));
    
    % How many of each neuron type are there?
    num.d1  = size(list.d1,  1);
    num.d2  = size(list.d2,  1);
    num.msn = size(list.msn, 1);
    num.fsi = size(list.fsi, 1);
    num.gap = size(connections.gap, 1);

    % (Approximate) number of MSNs of each type and FSIs to leave as background only
    num.bkg_msn = round(num.msn * (conn.bkg_msn / 100) / 2);
    num.bkg_fsi = round(num.fsi * (conn.bkg_fsi / 100));

    % Number of MSNs of each type and FSIs to put in each channel
    % num.msn_ch = floor((num.msn / 2 - num.bkg) / conn.ch_all);
    num.d1_ch  = num.d1 / conn.ch_all;
    num.d2_ch  = num.d1 / conn.ch_all;
    
    num.d1_in  = (num.d1 - num.bkg_msn) / conn.ch_all;
    num.d2_in  = (num.d2 - num.bkg_msn) / conn.ch_all;
    num.fsi_in = num.fsi - num.bkg_fsi;
    
    % TODO: SPLIT INPUT LIST GENERATION INTO TWO FUNCTIONS FOR PHYS / STAT?
        
    %% INPUT CONNECTIONS
    if flags.progress
        timer.conn1 = tic;
        fprintf('\nCreating connection lists:\n');
        fprintf('1) Cortical channel connections… ')       
    end
        
    if conn.ch_all == 1
        % Create connection ID for current background / channel width profile
        attr.conn_id = matlab.lang.makeValidName(strcat(conn_prefix, '_1CH'));

        % In the single-channel model it doesn't matter which neurons don't
        % receive connections since neuron ID is not associated with location
        output_lists.ch1.d1.(attr.conn_id)  = [ 0 : num.d1_in  - 1                          ; 0 : num.d1_in  - 1]';
        output_lists.ch1.d2.(attr.conn_id)  = [(0 : num.d2_in  - 1) + num.d1_in             ; 0 : num.d2_in  - 1]';
        output_lists.ch1.fsi.(attr.conn_id) = [(0 : num.fsi_in - 1) + num.d1_in + num.d2_in ; 0 : num.fsi_in - 1]';
        
    elseif conn.ch_all == 2       
        % Physically partition a two-channel striatum      
        for i = 1:2
            % Set dynamic structure fieldname
            msn = sprintf('d%d', i);


            % MAYBE CHECK FOR PHYS HERE
            % COULD DO CHANNEL WIDTH FOR STAT MODEL - HOW MANY NEURONS END UP IN BOTH CHANNELS
            if flags.phys_ch 
                % Create connection ID for current background / channel width profile
                attr.conn_id = matlab.lang.makeValidName(strcat(conn_prefix, '_wCH', num2str(phys.ch_width), '_2CH'));
                
                % MSNs assigned to channel 1 or 2 based on physical location on striatal X-axis, modified based on width value.
                % Width <0.5 creates a background-only gap between channels, width >0.5 creates a region with MSNs in both channels   
                temp.ch1.(msn) = list.(msn)(striatum.neurons(list.(msn)(:, 1), 1) <= 0         + (phys.size * (phys.ch_width / 100)), :);
                temp.ch2.(msn) = list.(msn)(striatum.neurons(list.(msn)(:, 1), 1) >= phys.size - (phys.size * (phys.ch_width / 100)), :);
            else 
                % Create connection ID for current background / channel width profile
                attr.conn_id = matlab.lang.makeValidName(strcat(conn_prefix, '_2CH'));
                
                % TODO: PROCEDURE FOR NON PARTITIONED 2CH LISTS
                % TODO: FIX THIS!!!
%                 ch = sprintf('d%d_ch', i);
                
                temp.ch1.(msn) = list.(msn)(1 : (num.(msn) / 2), :);
                temp.ch2.(msn) = list.(msn)((num.(msn) / 2) + 1 : num.(msn), :);           
            end

            for j = 1:conn.ch_all
                % Set dynamic structure fieldname
                ch = sprintf('ch%d', j);

                % Trim MSN and FSI lists according to requested background-only percentage
                temp.(ch).(msn) = temp.(ch).(msn)(1:end - (floor(size(temp.(ch).(msn), 1) * (conn.bkg_msn / 100))), :);
                temp.(ch).fsi   = list.fsi(       1:end - (floor(size(list.fsi, 1)        * (conn.bkg_fsi / 100))), :);

                % Channel connections to striatum
                % FROM:  Each cortical channel (-1 for SpineCreator 0-indexing)
                % TO:    Each MSN or FSI in each channel
                % DELAY: N/A                  
                if i==2
                    % D2 connections come after D1 connections
                    output_lists.(ch).(msn).(attr.conn_id) = ...
                        [(0 : length(temp.(ch).(msn)) - 1) + length(temp.(ch).d1) ; temp.(ch).(msn)(:,2)']';

                    % FSI connections come after all MSN connections
                    output_lists.(ch).fsi.(attr.conn_id) = ...
                        [(0 : length(temp.(ch).fsi) - 1) + length(temp.(ch).d1) + length(temp.(ch).d2) ; temp.(ch).fsi(:,2)']';
                else
                    % No special considerations for D1 connections
                    output_lists.(ch).(msn).(attr.conn_id) = ...
                        [0 : length(temp.(ch).(msn)) - 1 ; temp.(ch).(msn)(:,2)']';
                end                      
            end
        end
    else
        error('Three or more channels currently not supported')
    end
    
    % Save connection lists
    if flags.save
        if ~exist(fullfile(attr.list_path, attr.conn_id), 'dir')
            mkdir(fullfile(attr.list_path, attr.conn_id))
        end
        
        for i = 1:conn.ch_all
            % Set dynamic structure fieldname
            ch = sprintf('ch%d', i);       
            name.src = sprintf('CH%d_input', i);

            % For both D1 and D2 MSNs
            for j = 1:2           
                d_dst = sprintf('d%d', j);
                name.dst = sprintf('Striatum_D%d', j);

                % Create both AMPA and NMDA connections to MSNs
                for k = 0:1  
                    name.syn = sprintf('syn%d', k); 
                    save_list(fullfile(attr.list_path, attr.conn_id), output_lists.(ch).(d_dst).(attr.conn_id), name, flags);
                end
            end

            % FSIs only use a single synapse
            name.dst = 'Striatum_FSI';
            name.syn = 'syn0';
            save_list(fullfile(attr.list_path, attr.conn_id), output_lists.(ch).fsi.(attr.conn_id), name, flags);
        end
    end
    
    if flags.progress
    	fprintf('done! (%1.2fs)\n', toc(timer.conn1))
    end
    
    %% INTRASTRIATAL CONNECTIONS
    if flags.progress
        timer.conn2 = tic;
        fprintf('2) Striatal connections… ')
    end
        
    % From D1 / D2 MSNs
    for i = 1:2
        msn_i = sprintf('d%d', i); 
        name.src = sprintf('Striatum_D%d', i);

        % To D1 / D2 MSNs
        for j = 1:2 
            msn_j = sprintf('d%d',j);
            name.dst = sprintf('Striatum_D%d', j);
            name.syn = 'syn0';
            
            % If connections already exist, don't recreate
            if ~exist(fullfile(attr.list_path, ['conn_', name.src, '_to_', name.dst, '_syn0.csv']), 'file')        
                output_lists.(msn_i).(msn_j) = create_list(connections.(msn_i).(msn_j), list.(msn_i), list.(msn_j), phys.cv_msnmsn);
                
                % Save connection lists
                if flags.save
                    save_list(attr.list_path, output_lists.(msn_i).(msn_j), name, flags);
                end
            end
        end
                 
        % From FSIs to MSNs
        name.src = 'Striatum_FSI';
        name.dst = sprintf('Striatum_D%d', i);
        name.syn = 'syn0';
        
        % If connections already exist, don't recreate
        if ~exist(fullfile(attr.list_path, ['conn_', name.src, '_to_', name.dst, '_syn0.csv']), 'file')            
            output_lists.fsi.(msn_i) = create_list(connections.fsi.(msn_i), list.fsi, list.(msn_i), phys.cv_fsimsn);
            
            % Save connection lists
            if flags.save
                save_list(attr.list_path, output_lists.fsi.(msn_i), name, flags);
            end
        end
    end

    % From FSIs to FSIs (GABA)
    name.src = 'Striatum_FSI';
    name.dst = 'Striatum_FSI';
    name.syn = 'syn0';

    % If connections already exist, don't recreate
    if ~exist(fullfile(attr.list_path, ['conn_', name.src, '_to_', name.dst, '_syn0.csv']), 'file')        
        output_lists.fsi.fsi = create_list(connections.fsifsi, list.fsi, list.fsi, phys.cv_fsifsi);
        
        % Save connection lists
        if flags.save
            save_list(attr.list_path, output_lists.fsi.fsi, name, flags);
        end

        % From FSIs to FSIs (Gap)
        % Convert MatLab neuron IDs to SpineCreator IDs
        try
            [~, src] = ismember(connections.gap(:, 1), list.fsi(:, 1));
            [~, dst] = ismember(connections.gap(:, 2), list.fsi(:, 1));
        catch 
        end

        % Create list of FSI-FSI gap connections
        % Gap junctions have a nonstandard format
        try
            output_lists.gap.in1 =  [list.fsi(src, 2)' ; 0 : num.gap - 1  ]';
            output_lists.gap.in2 =  [list.fsi(dst, 2)' ; 0 : num.gap - 1  ]';
            output_lists.gap.out1 = [0 : num.gap - 1   ; list.fsi(src, 2)']';
            output_lists.gap.out2 = [0 : num.gap - 1   ; list.fsi(dst, 2)']';        
        catch
            error('Could not create gap junctions!')
        end
        
        % Save connection lists
        if flags.save
            for i = 0:1
                g_in = sprintf('in%d', i + 1); 
                g_out = sprintf('out%d', i + 1); 

                name.syn = sprintf('syn%d', i);

                name.src = 'Striatum_FSI';               
                name.dst = 'FSI_GAP';
                save_list(attr.list_path, output_lists.gap.(g_in), name, flags);

                name.src = 'FSI_GAP';
                name.dst = 'Striatum_FSI';
                save_list(attr.list_path, output_lists.gap.(g_out), name, flags);        
            end
        end
    end
    
    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn2))
    end
    
    %% FUNCTIONS
    function conn_list = create_list(input_list, src_list, dst_list, varargin)
        % Convert MatLab neuron IDs to SpineCreator IDs and create connection list
        [~, s] = ismember(input_list(:, 1), src_list(:, 1));
        [~, d] = ismember(input_list(:, 2), dst_list(:, 1));

%         % If a connection delay multiplier exists, use it
%         try
%             conn_list = [src_list(s, 2)' ; dst_list(d, 2)' ; input_list(:, 3)' .* attr.conn.delay_mult]'; 
% %             conn_list = [src_list(s, 2)' ; dst_list(d, 2)' ; input_list(:, 3)' .* 0.02 .* 0.95 + 0.1]'; 
%             
%         catch
%             conn_list = [src_list(s, 2)' ; dst_list(d, 2)' ; input_list(:, 3)']';
%         end

%         if size(input_list, 2) == 2
%             % No third column == no distance == statistical striatum, use standard delay for all connections
%             conn_list = [src_list(s, 2)' ; dst_list(d, 2)' ; repmat(attr.conn.delay_min * attr.conn.delay_mult, 1, size(input_list, 1))]';
%         elseif size(input_list, 2) == 3
%             % Third column present == physical distance in μm
%             
%             conn_list = [src_list(s, 2)' ; dst_list(d, 2)' ; max(input_list(:, 3) ./ (phys.cv_msnmsn * 1000 / attr.conn.delay_mult), conn.delay_min)']';
%         else
%             error('Incorrectly sized connection list')               
%         end
        
        % If a connection distance and conductance velocity are given, construct
        % distance-based connection delays
        if size(input_list, 2) == 3 && exist('varargin', 'var')
            velocity = varargin{1};
            
            % Multiply conduction velocity in m/s by 1,000 for equivalent in μm/ms
            conn_list = [src_list(s, 2)' ; dst_list(d, 2)' ; max(input_list(:, 3) ./ (velocity * 1000 / conn.delay_mult), conn.delay_min)']';
        else
            conn_list = [src_list(s, 2)' ; dst_list(d, 2)' ; repmat(conn.delay_min * conn.delay_mult, 1, size(input_list, 1))]';
        end
            
    end  
end


%% FUNCTIONS
function[] = save_list(path, file, name, flags)
    % Given a pathname, a connection list and name, this will save the
    % connection list to a CSV and optionally convert it to binary format for
    % direct import to SpineCreator

    % Save connection list to CSV
    filename = ['conn_', name.src, '_to_', name.dst, '_', name.syn];
    fid = fopen(fullfile(path, [filename, '.csv']), 'w');

    % With or without connection delays as appropriate
    if size(file,2) == 3
        fprintf(fid, '%d, %d, %1.1f\r\n', transpose(file));
    elseif size(file,2) == 2
        fprintf(fid, '%d, %d\r\n', transpose(file));
    else
        error('Incorrect number of columns in connection %s', filename);
    end
    fclose(fid);

    % Export connection list as binary file if needed
    if flags.binary
        bin_convert(file, fullfile(path, filename));
    end


    function[] = bin_convert(data, filename)
        % Given a two- or three-column array of data representing input neuron,
        % output neuron and (optionally) delay, this function converts the data to 
        % binary format so that it can be dropped straight into the SpineCreator model directory.

        intArray = int32([data(1 : end, 1)' ; data(1 : end, 2)']');
        if size(data, 2) == 3
            floatArray = single(data(1 : end, 3));
        end

        fileID = fopen([filename, '.bin'], 'w');

        if exist('floatArray', 'var')
            for n = 1:size(data, 1)
                fwrite(fileID, intArray(n, :), 'int32');
                fwrite(fileID, floatArray(n), 'single');
            end
        else
            for n = 1:size(data,1)
                fwrite(fileID, intArray(n, :), 'int32');
            end
        end

        fclose(fileID);
    end       
end