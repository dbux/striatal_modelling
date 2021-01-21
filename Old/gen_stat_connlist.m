function[connections] = gen_stat_connlist(attr, conn, stat, flags)
    % Takes attributes and flags and creates statistical striatal connection lists

    % Start list timer
    timer.list = tic;

    % Ensure we have a round number of MSNs
    if mod(stat.num_msn, 2) ~= 0
        stat.num_msn = stat.num_msn + 1;
        fprintf('Odd number of MSNs detected! Have a bonus one!\n');
    end

    % Number of each type of MSN
    num_d = stat.num_msn / 2;         

    % Number of background-only MSNs
    num_bkg = stat.num_msn * (conn.bkg_msn / 100);       

    % Number of MSNs of each type in each channel
    num_ch = floor(((stat.num_msn - num_bkg) / (conn.ch_all * 2)));

    % Given the number of neurons and the number of expected connections, what
    % is the probability of a connection between any two neurons?
    prob_msnmsn = stat.con_msnmsn / stat.num_msn;
    prob_fsimsn = stat.con_fsimsn / stat.num_msn;
    prob_fsifsi = stat.con_fsifsi / stat.num_fsi;
    prob_fsigap = stat.con_fsigap / stat.num_fsi;

    % Directory for saving connection lists
    strpath = ...
        [attr.path num2str(datestr(now, 'yy.mm.dd_HH.MM')) '_' ...
        num2str(stat.num_msn) '+' num2str(stat.num_fsi) '_' num2str(conn.ch_all) 'CH_' num2str(conn.bkg_msn) 'BKG_stat/'];
    listpath = [strpath 'connection_lists/'];
    mkdir(strpath);
    mkdir(listpath);

    %% CORTICAL CHANNEL CONNECTIONS
    if flags.progress
        fprintf('\nCreating statistical connection lists:\n');
        fprintf('1) Cortical channel connections... ')
    end
    timer.conn1 = tic;

    % Position of the first and last neuron in the channel
    ch_first = 0;
    ch_last = num_ch - 1;
    ch_next = ch_last + num_ch;

    for i = 1:conn.ch_all
        % Channel variable names
        ch = sprintf('ch%d', i);       
        name.src = sprintf('CH%d_input', i);

        %Cortical inputs to FSI and BG use synapse 0
        name.syn = 'syn0';

        %% CORTEX-FSI CONNECTIONS
        connections.cortex.(ch).fsi = [0 : stat.num_fsi - 1 ; 0 : stat.num_fsi - 1]';

        %% CORTEX-BG CONNECTIONS
        if num_ch > 500
            connections.cortex.(ch).bg = [0 : 499 ; zeros(1, 500) + (i - 1)]';
        else
            connections.cortex.(ch).bg = [0 : num_ch - 1 ; zeros(1, num_ch) + (i - 1)]';
        end

        % Save connection lists
        if flags.save
            % Cortex-FSI
            name.dst = 'Striatum_FSI'; 
            save_list(listpath, connections.cortex.(ch).fsi, name, flags);

            % Cortex-MCtx
            name.dst = 'Motor_Cortex';
            save_list(listpath, connections.cortex.(ch).bg, name, flags);

            % Cortex-STN
            name.dst = 'STN';
            save_list(listpath, connections.cortex.(ch).bg, name, flags); 
        end

        %% CORTEX-STRIATUM CONNECTIONS
        for j = 1:2
            % Channel variable names
            cor_d = sprintf('d%d', j);           
            name.dst = sprintf('Striatum_D%d', j);

            % Create one-to-one connection between cortical inputs and MSNs in channel
            connections.cortex.(ch).(cor_d) = [0 : num_ch - 1 ; ch_first : ch_last]';

            % Save connection lists
            if flags.save
                % Create both AMPA and NMDA connections to MSNs
                for k = 0:1
                    name.syn = sprintf('syn%d', k);                
                    save_list(listpath, connections.cortex.(ch).(cor_d), name, flags);
                end
            end
        end    

        % Increment channel start/end markers for the next channel
        ch_first = ch_last + 1;
        ch_last = ch_next;
        ch_next = ch_last + num_ch;     
    end

    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn1))
    end

    %% CORTICAL BACKGROUND CONNECTIONS
    if flags.progress
        fprintf('2) Cortical background connections... ')
    end
    timer.conn2 = tic;

    if conn.bkg_msn > 0  
        %% CORTEX-BG BACKGROUND CONNECTIONS
        connections.cortex.bkg.bg = [];
        for i = 1:attr.ch_all
            if num_ch > attr.max_bg
                connections.cortex.bkg.bg = [connections.cortex.bkg.bg ; [0 : attr.max_bg - 1 ; zeros(1, attr.max_bg) + (i - 1)]'];
            else
                connections.cortex.bkg.bg = [connections.cortex.bkg.bg ; [0 : num_ch - 1      ; zeros(1, num_ch) + (i - 1)]'     ];
            end
        end

        %% CORTEX-STRIATUM BACKGROUND CONNECTIONS
        connections.cortex.bkg.d1 = [0 : num_d - 1        ; 0 : num_d - 1       ]';
        connections.cortex.bkg.d2 = [0 : num_d - 1        ; 0 : num_d - 1       ]';
        connections.cortex.bkg.fsi= [0 : stat.num_fsi - 1 ; 0 : stat.num_fsi - 1]';

        % Save connection lists
        if flags.save   
            name.src = 'BKG_input';

            % For both D1 and D2 MSNs
            for i = 1:2
                d_dst = sprintf('d%d', i);
                name.dst = sprintf('Striatum_D%d', i);

                % Create both AMPA and NMDA connections to MSNs
                for j = 0:1
                    name.syn = sprintf('syn%d', j);               
                    save_list(listpath, connections.cortex.bkg.(d_dst), name, flags);
                end
            end

            % Other populations only use a single synapse
            name.syn = 'syn0';

            % Cortex-FSI
            name.dst = 'Striatum_FSI';
            save_list(listpath, connections.cortex.bkg.fsi, name, flags);

            % Cortex-MCtx
            name.dst = 'Motor_Cortex';
            save_list(listpath, connections.cortex.bkg.bg, name, flags);

            % Cortex-STN
            name.dst = 'STN';
            save_list(listpath, connections.cortex.bkg.bg, name, flags);
        end
    end

    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn2))
    end

    %% STRIATAL GABA CONNECTIONS
    if flags.progress
        fprintf('3) Striatal GABA connections... ')
    end
    timer.conn3 = tic;

    % GABA projections always use syn0
    name.syn = 'syn0';
    for i = 1:2       
        %% MSN-MSN GABA CONNECTIONS
        % Source MSN type
        con_s = sprintf('d%d', i); 
        name.src = sprintf('Striatum_D%d', i);
        for j = 1:2        
            % Destination MSN type
            con_d = sprintf('d%d', j);
            name.dst = sprintf('Striatum_D%d', j);

            % Create an all-to-all connection list between D(i) and D(j) MSNs
            connections.striatum.(con_s).(con_d) = gen_list_allall(num_d, num_d);

            % Trim the connections to match the expected probability value
            connections.striatum.(con_s).(con_d) = ...
                connections.striatum.(con_s).(con_d)((connections.striatum.(con_s).(con_d)(:,3) <= prob_msnmsn),:);

            % Replace the random probability value with a fixed connection delay
    %         connections.striatum.(con_s).(con_d)(:,3) = ...
    %             rand(length(connections.striatum.(con_s).(con_d)), 1) * attr.delay_mult + conn.delay_min;
             connections.striatum.(con_s).(con_d)(:,3) = conn.delay_min;

            % Save connection lists
            if flags.save
                save_list(listpath, connections.striatum.(con_s).(con_d), name, flags);
            end
        end

        %% FSI-MSN CONNECTIONS   
        name.src = 'Striatum_FSI';
        name.dst = sprintf('Striatum_D%d', i);

        % Create all-to-all FSI-MSN connection list
        connections.striatum.fsi.(con_s) = gen_list_allall(stat.num_fsi, num_d);

        % If FSI-MSN probability is greater than one (i.e. more than one
        % connection should be made) then connections with a probability value
        % less than (probability - 1) are duplicated and none are removed
        if prob_fsimsn > 1
            prob_curr = prob_fsimsn;
            while prob_curr > 1
                fsi_extra = connections.striatum.fsi.(con_s)((connections.striatum.fsi.(con_s)(:,3) <= prob_curr - 1),:);
                connections.striatum.fsi.(con_s) = [connections.striatum.fsi.(con_s) ; fsi_extra];
                prob_curr = prob_curr - 1;
            end
        else
            % Trim the connections to match the expected probability value
            connections.striatum.fsi.(con_s) = connections.striatum.fsi.(con_s)((connections.striatum.fsi.(con_s)(:,3) <= prob_fsimsn),:);
        end

        % Replace the random probability value with a fixed connection delay
    %     connections.striatum.fsi.(con_s)(:,3) = rand(length(connections.striatum.fsi.(con_s)), 1) * attr.delay_mult + conn.delay_min;
        connections.striatum.fsi.(con_s)(:,3) = conn.delay_min;

        % Save connection lists
        if flags.save
            save_list(listpath, connections.striatum.fsi.(con_s), name, flags);
        end
    end

    %% FSI-FSI CONNECTIONS
    name.src = 'Striatum_FSI';
    name.dst = 'Striatum_FSI';

    % Create all-to-all FSI-FSI connection list
    connections.striatum.fsi.fsi = gen_list_allall(stat.num_fsi, stat.num_fsi);

    % Trim the connections to match the expected probability value
    connections.striatum.fsi.fsi = connections.striatum.fsi.fsi((connections.striatum.fsi.fsi(:,3) <= prob_fsifsi),:);

    % Replace the random probability value with a fixed connection delay
    % connections.striatum.fsi.fsi(:,3) = rand(size((connections.striatum.fsi.fsi), 1), 1) * attr.delay_mult + conn.delay_min;
    connections.striatum.fsi.fsi(:,3) = conn.delay_min;

    % Save connection lists
    if flags.save
        save_list(listpath, connections.striatum.fsi.fsi, name, flags);
    end

    %% FSI GAP CONNECTIONS
    % The procedure to create FSI gap connections is a bit different:
    % Gap connections can only exist where a FSI-FSI connection is already present
    connections.striatum.fsi.gap = connections.striatum.fsi.fsi;

    % Replace the random connection delay with a random probability value
    connections.striatum.fsi.gap(:,3) = rand(1, size(connections.striatum.fsi.fsi, 1)');

    % Trim the connections to match the expected probability value
    connections.striatum.fsi.gap = connections.striatum.fsi.gap((connections.striatum.fsi.gap(:,3) <= prob_fsigap),:);

    % Gap connections have a connection delay of 0
    % connections.striatum.fsi.gap(:,3) = 0;
    connections.striatum.fsi.gap(:,3) = [];

    % Because gap connections have their own 'neuron' population each connection
    % must go via an intermediate, zero-indexed list
    num_gap = size(connections.striatum.fsi.gap, 1);

    connections.striatum.fsi.gap_in1 = [connections.striatum.fsi.gap(:,1)' ; 0:num_gap - 1']';
    connections.striatum.fsi.gap_in2 = [connections.striatum.fsi.gap(:,2)' ; 0:num_gap - 1']';

    connections.striatum.fsi.gap_out1 = [0:num_gap - 1' ; connections.striatum.fsi.gap(:,1)']';
    connections.striatum.fsi.gap_out2 = [0:num_gap - 1' ; connections.striatum.fsi.gap(:,2)']';

    % Save connection lists
    if flags.save
        for i = 0:1
            g_in = sprintf('gap_in%d', i + 1); 
            g_out = sprintf('gap_out%d', i + 1); 

            name.syn = sprintf('syn%d', i);

            name.src = 'Striatum_FSI';               
            name.dst = 'FSI_GAP';
            save_list(listpath, connections.striatum.fsi.(g_in), name, flags);

            name.src = 'FSI_GAP';
            name.dst = 'Striatum_FSI';
            save_list(listpath, connections.striatum.fsi.(g_out), name, flags);        
        end
    end

    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn3))
    end

    %% STRIATAL PEPTIDE CONNECTIONS
    if flags.progress
        fprintf('4) Striatal peptide connections... ')
    end
    timer.conn4 = tic;

    % Peptide projections always use syn1
    name.syn = 'syn1';
    for i = 1:2
        % Source MSN type
        d_src = sprintf('d%d', i);
        name.src = sprintf('Striatum_D%d', i);

        for j = 1:2 
            % Destination MSN type
            d_dst = sprintf('d%d',j);
            name.dst = sprintf('Striatum_D%d', j);

            % Empty arrays for diffuse, unidirectional and pruned connection lists
            connections.dff.(d_src).(d_dst) = [];
            connections.uni.(d_src).(d_dst) = [];       
            connections.prn.(d_src).(d_dst) = [];

            % Position of the first and last neuron in the channel
            seq_ch_first = 0;
            seq_ch_last = num_ch - 1;
            seq_ch_next = seq_ch_last + num_ch;

            % For each action channel
            for k = 1:attr.ch_seq
                ch_src = sprintf('ch%d', k);

                % Array of neurons in current channel
                seq_mem = seq_ch_first:seq_ch_last;

                % Position of the first and last neuron in the channel
                ch_first = 0;
                ch_last = num_ch - 1;
                ch_next = ch_last + num_ch;

                for l = 1:attr.ch_all
                    ch_dst = sprintf('ch%d', l);

                    % Array of neurons in current channel
                    ch_mem = ch_first:ch_last;

                    connections.(d_src).(ch_src).(d_dst).(ch_dst) = ...
                        connections.striatum.(d_src).(d_dst)...
                        (ismember(connections.striatum.(d_src).(d_dst)(:,1), seq_mem) & ...
                        ismember(connections.striatum.(d_src).(d_dst)(:,2), ch_mem), :);

                    % Diffuse and pruned connections include all connections
                    connections.dff.(d_src).(d_dst) = [connections.dff.(d_src).(d_dst) ; connections.(d_src).(ch_src).(d_dst).(ch_dst)];
                    connections.prn.(d_src).(d_dst) = [connections.prn.(d_src).(d_dst) ; connections.(d_src).(ch_src).(d_dst).(ch_dst)];

                    % Unidirectional connection if source channel < number of
                    % channels in sequence AND dest. channel = source channel + 1
                    if l == (k+1) && k < attr.ch_seq
                        connections.uni.(d_src).(d_dst) = [connections.uni.(d_src).(d_dst) ; connections.(d_src).(ch_src).(d_dst).(ch_dst)];
                    end  

                    % Connections to be pruned if source MSN is D1 and prune
                    % source and destination channels match those set by user
                    if i == 1 && k == attr.prn_src && l == attr.prn_dst                                      
                        connections.prn.(d_src).(d_dst)...
                            (ismember(connections.prn.(d_src).(d_dst), connections.(d_src).(ch_src).(d_dst).(ch_dst), 'rows'),:) = [];
                    end

                    % Increment channel start/end markers for the next channel
                    ch_first = ch_last + 1;
                    ch_last = ch_next;
                    ch_next = ch_last + num_ch; 
                end

                % Increment channel start/end markers for the next channel
                seq_ch_first = seq_ch_last + 1;
                seq_ch_last = seq_ch_next;
                seq_ch_next = seq_ch_last + num_ch;
            end

            % Save connection lists
            if flags.save
                if i == 1
                    name.dst = sprintf('Striatum_D%d (Diffuse)', j);
                    save_list(listpath, connections.dff.(d_src).(d_dst), name, flags);

                    name.dst = sprintf('Striatum_D%d (Pruned)', j);
                    save_list(listpath, connections.prn.(d_src).(d_dst), name, flags);

                    name.dst = sprintf('Striatum_D%d (Unidirectional)', j);
                    save_list(listpath, connections.uni.(d_src).(d_dst), name, flags);
                else
                    name.dst = sprintf('Striatum_D%d', j);
                    save_list(listpath, connections.dff.(d_src).(d_dst), name, flags);
                end
            end
        end
    end              

    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn4))
    end

    %% STRIATUM-GP (& MCtx-MSN) CONNECTIONS
    if flags.progress
        fprintf('5) Striatal GP connections... ')
    end
    timer.conn5 = tic;

    for i = 1:2 
        con_s = sprintf('d%d', i);

        % Position of the first and last neuron in the channel
        ch_first = 0;
        ch_last = num_ch - 1;
        ch_next = ch_last + num_ch;

        % From striatum to BG loop and from MCtx to striatum   
        connections.striatum.(con_s).bg = [];
        connections.striatum.bg.(con_s) = [];
        connections.striatum.bg.fsi = [];

        for j = 1:attr.ch_all          
            con_d = sprintf('gp%d', j); 

            % Limit striatum-GP connections if necessary
            if num_ch > attr.max_bg
                connections.striatum.(con_s).bg = ...
                    [connections.striatum.(con_s).bg ; [ch_first : (ch_first + (attr.max_bg - 1)) ; zeros(1, attr.max_bg) + (j - 1)]'];
            else
                connections.striatum.(con_s).bg = ...
                    [connections.striatum.(con_s).bg ; [ch_first : ch_last                        ; zeros(1, num_ch) + (j - 1)]'     ];
            end

            connections.striatum.bg.(con_s) = [connections.striatum.bg.(con_s) ; [zeros(1, num_ch) + (j - 1)       ; ch_first : ch_last  ]'];
            connections.striatum.bg.fsi     = [connections.striatum.bg.fsi     ; [zeros(1, stat.num_fsi) + (j - 1) ; 0 : stat.num_fsi - 1]'];     

            % Increment channel start/end markers for the next channel
            ch_first = ch_last + 1;
            ch_last = ch_next;
            ch_next = ch_last + num_ch;
        end

        % Save connection lists
        if flags.save
            switch i
                case 1
                    % D1 MSN to GPi/SNr
                    name.dst = 'GPi_SNr';
                case 2
                    % D2 MSN to GPe
                    name.dst = 'GPe';
            end   

            % Striatum to BG
            name.src = sprintf('Striatum_D%d', i);
            name.syn = 'syn0';  
            save_list(listpath, connections.striatum.(con_s).bg, name, flags);

            % MCtx to FSI
            name.src = 'MCtx_R2S';
            name.dst = 'Striatum_FSI';
            name.syn = 'syn0';
            save_list(listpath, connections.striatum.bg.fsi, name, flags);

            % MCtx to MSNs
            name.dst = sprintf('Striatum_D%d', i);
            for j = 0:1
                name.syn = sprintf('syn%d', j);
                save_list(listpath, connections.striatum.bg.(con_s), name, flags); 
            end
        end
    end

    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn5))
    end


    if flags.progress
        fprintf('All connection lists generated in %1.2f minutes.\n', toc(timer.list) / 60)
    end

    % Save connections file
    filename = [strpath 'connections.mat'];
    save(filename, 'connections', '-v7.3');

    %% CONNECTION STATISTICS
    fprintf('\n* Connection lists generated and saved *\n');
    fprintf('\n- General statistics -\n');
    fprintf('Number of MSNs: %d\n', stat.num_msn);
    fprintf('Number of FSIs: %d\n', stat.num_fsi);
    fprintf('Number of action channels: %d (including sequence of %d)\n', attr.ch_all, attr.ch_seq);

    fprintf('\n- MSN connection statistics -\n');
    % Get MSNs - 1 MSN
    msn_count = 0;
    for i = 1:2
        for j = 1:2
            con_s = sprintf('d%d', i);
            con_d = sprintf('d%d', j); 
            % Adapted from http://uk.mathworks.com/matlabcentral/answers/19042-finding-duplicate-values-per-column#answer_25447
            uniqueX = unique(connections.striatum.(con_s).(con_d)(:,2));
            countOfX = hist(connections.striatum.(con_s).(con_d)(:,2),uniqueX);
            indexToRepeatedValue = (countOfX~=1);
            numberOfAppearancesOfRepeatedValues = countOfX(indexToRepeatedValue);
            msn_count = msn_count + mean(numberOfAppearancesOfRepeatedValues);
        end
    end
    fprintf('D1-D1: %d\n', length(connections.striatum.d1.d1)); 
    fprintf('D1-D2: %d\n', length(connections.striatum.d1.d2)); 
    fprintf('D2-D1: %d\n', length(connections.striatum.d2.d1)); 
    fprintf('D2-D2: %d\n', length(connections.striatum.d2.d2));
    fprintf('Expected (for any MSN-MSN group): %.1f\n', prob_msnmsn * num_d^2);
    fprintf('%1.3f MSNs-1 MSN (Expected: %d)\n', msn_count / 2, stat.con_msnmsn);

    fprintf('\n- FSI connection statistics -\n');
    % Get 1 FSI - MSNs
    fsim_count(1) = 0;
    fsim_count(2) = 0;
    for i = 1:2
        con_d = sprintf('d%d', i); 
        % Adapted from http://uk.mathworks.com/matlabcentral/answers/19042-finding-duplicate-values-per-column#answer_25447
        uniqueM = unique(connections.striatum.fsi.(con_d)(:,1));
        uniqueF = unique(connections.striatum.fsi.(con_d)(:,2));

        countOfM = hist(connections.striatum.fsi.(con_d)(:,1),uniqueM);
        countOfF = hist(connections.striatum.fsi.(con_d)(:,2),uniqueF);

        indexToRepeatedValueM = (countOfM~=1);
        indexToRepeatedValueF = (countOfF~=1);

        numberOfAppearancesOfRepeatedValuesM = countOfM(indexToRepeatedValueM);
        numberOfAppearancesOfRepeatedValuesF = countOfF(indexToRepeatedValueF);

        fsim_count(1) = fsim_count(1) + mean(numberOfAppearancesOfRepeatedValuesM);
        fsim_count(2) = fsim_count(2) + mean(numberOfAppearancesOfRepeatedValuesF);
    end
    fsim_count(2) = fsim_count(2) / 2;

    fprintf('FSI-D1: %d\n', length(connections.striatum.fsi.d1)); 
    fprintf('FSI-D2: %d\n', length(connections.striatum.fsi.d2)); 
    fprintf('Expected (for any FSI-MSN group): %.1f\n', prob_fsimsn * stat.num_fsi * num_d);
    fprintf('1 FSI-%1.3f MSNs (Expected: %d)\n', fsim_count(1), stat.con_fsimsn);
    fprintf('%1.3f FSIs-1 MSN (Expected: 30.6)\n', fsim_count(2));

    % Get FSIs - 1 FSI
    % Adapted from http://uk.mathworks.com/matlabcentral/answers/19042-finding-duplicate-values-per-column#answer_25447
    uniqueX = unique(connections.striatum.fsi.fsi(:,1));
    countOfX = hist(connections.striatum.fsi.fsi(:,1),uniqueX);
    indexToRepeatedValue = (countOfX~=1);
    numberOfAppearancesOfRepeatedValues = countOfX(indexToRepeatedValue);
    fsif_count = mean(numberOfAppearancesOfRepeatedValues);
    fprintf('\nFSI-FSI: %d\n', length(connections.striatum.fsi.fsi)); 
    fprintf('Expected: %.1f\n', prob_fsifsi * stat.num_fsi^2);
    fprintf('%1.3f FSIs-1 FSI (Expected: %1.1f)\n', fsif_count, stat.con_fsifsi);

    fprintf('\nFSI gap: %d\n', length(connections.striatum.fsi.gap)); 
    fprintf('Expected: %.1f\n', prob_fsigap * length(connections.striatum.fsi.fsi));

    fprintf('\n- Neuropeptide statistics -\n');
    fprintf('Substance P: %d%% unidirectional bias\n', attr.uni_bias);
    % fprintf('Unidirectional connections: %d to D1, %d to D2\n', ...
    %     sum(ismember([sp.d1.final(:,1), sp.d1.final(:,2)], [sp.d1.uni(:,1), sp.d1.uni(:,2)], 'rows') ~= 0),...
    %     sum(ismember([sp.d2.final(:,1), sp.d2.final(:,2)], [sp.d2.uni(:,1), sp.d2.uni(:,2)], 'rows') ~= 0));
    fprintf('Unidirectional connections: %d to D1, %d to D2\n', ...
        length(connections.uni.d1.d1), length(connections.uni.d1.d2));
    % fprintf('Diffuse connections: %d to D1, %d to D2\n', ...
    %     sum(ismember([sp.d1.final(:,1), sp.d1.final(:,2)], [sp.d1.diff(:,1), sp.d1.diff(:,2)], 'rows') ~= 0),...
    %     sum(ismember([sp.d2.final(:,1), sp.d2.final(:,2)], [sp.d2.diff(:,1), sp.d2.diff(:,2)], 'rows') ~= 0));
    fprintf('Diffuse connections: %d to D1, %d to D2\n', ...
        length(connections.dff.d1.d1), length(connections.dff.d1.d2));
end