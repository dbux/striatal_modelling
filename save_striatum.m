function save_striatum(striatum, connections, list, attr) 
% Export striatal co-ordinate data to a Veusz-compatible CSV

    % Extract necessary attributes
    conn = attr.conn;
    flags = attr.flags;
    
    if flags.progress
        timer.save = tic;
        fprintf('\nSaving neuron and connection data… ')
    end
    
    % Get physical co-ordinates of all MSNs and FSIs
    striatum.msn.d1  = striatum.neurons(list.d1(:, 1), :);
    striatum.msn.d2  = striatum.neurons(list.d2(:, 1), :);
    striatum.msn.all = [striatum.msn.d1 ; striatum.msn.d2];
    striatum.fsi.all = striatum.neurons(list.fsi(:, 1), :);

    striatum.fsi.active    = [];
    striatum.fsi.inactive  = [];
    striatum.msn.active.d1 = [];
    striatum.msn.active.d2 = [];

    for i = 1:conn.ch_all
        % Set dynamic structure fieldname
        ch = sprintf('ch%d', i);

        % Get active FSIs
        [~, idf] = ismember(connections.(ch).fsi.(attr.conn_id)(:, 2), list.fsi(:, 2));
        striatum.fsi.active = unique([striatum.fsi.active ; striatum.fsi.all(idf, :)], 'rows');

        striatum.msn.(ch).all = [];

        % Get active MSNs
        for j = 1:2   
            % Set dynamic structure fieldname
            dx = sprintf('d%d', j);

            [~, idm] = ismember(connections.(ch).(dx).(attr.conn_id)(:, 2), list.(dx)(:, 2));
            striatum.msn.(ch).(dx)   = striatum.neurons(list.(dx)(idm, 1), :);
            striatum.msn.(ch).all    = [striatum.msn.(ch).all    ; striatum.msn.(ch).(dx)];
            striatum.msn.active.(dx) = [striatum.msn.active.(dx) ; striatum.msn.(ch).(dx)];
        end
    end

    % Get inactive neurons
    striatum.msn.inactive.d1  = striatum.msn.d1(~ismember(striatum.msn.d1, striatum.msn.active.d1, 'rows'), :);
    striatum.msn.inactive.d2  = striatum.msn.d2(~ismember(striatum.msn.d2, striatum.msn.active.d2, 'rows'), :);
    striatum.msn.inactive.all = [striatum.msn.inactive.d1 ; striatum.msn.inactive.d2];
    striatum.fsi.inactive     = striatum.fsi.all(~ismember(striatum.fsi.all, striatum.fsi.active, 'rows'), :);

    % Export co-ordinates of all MSNs and FSIs
    save_dir = fullfile(attr.save_path, 'neuron_data', attr.conn_id);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    export_striatum(striatum.msn, 'striatum.msn', save_dir)
    export_striatum(striatum.fsi, 'striatum.fsi', save_dir)

    % Save list of neurons
    list_name = fullfile(save_dir, 'list.mat');   
    save(list_name, 'list');
    
    if flags.progress
        fprintf('took %1.2f minutes. All done!\n', toc(timer.save)/60)
    end
    
    %% FUNCTIONS
    function export_striatum(s_data, s_name, s_dir)
        fields = fieldnames(s_data);
        for idx = 1:numel(fields)           
            % Get full name of current structure location
            curr_name = [s_name, '.', fields{idx}];

            % Get contents of next structure field
            s_field = s_data.(fields{idx});
            if isstruct(s_field)
                % Recurse function until nested structures end
                export_striatum(s_field, curr_name, s_dir);
            else
                % Export current co-ordinate matrix to disk
                cell2csv([fullfile(s_dir, curr_name), '.csv'], ...
                    [{[curr_name, '_X'], [curr_name, '_Y'], [curr_name, '_Z']} ; num2cell(s_field)]);
            end
        end   
    end
end