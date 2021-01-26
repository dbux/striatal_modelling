% Generates a 'physical' striatal microcircuit and associated connection 
% lists based on Humphries, Wood & Gurney (2010)

%% TODO LIST
% Move functions into main code section if only used once
% Consider speeding up generation by placing all neuons at once and performing distance checks later
% Add more detailed header information to 3D neuron data export

% Stop saving separate list.mat

%% PREAMBLE
% Reset initial state
clear variables; clc
% Specify RNG seed
rng(849)

% Add MATLAB path for HPC execution
addpath(genpath('~/MatLab'));
addpath(genpath('/home/ac1drb/MatLab'));

% timer.all = tic;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% Optionally define an existing striatum to generate connection lists for
s_ID = '21.01.25_20.06_84900+845';
% s_ID = '21.01.26_09.51_6000+60';

% Import attributes 
get_attributes


%% START
if attr.flags.progress
    timer.all = tic;
end
            
if ~exist('s_ID', 'var')
    if attr.flags.physical            
        % CREATE PHYSICAL STRIATUM
        striatum = gen_phys_striatum(attr);
         
        % GENERATE PHYSICAL INTRASTRIATAL CONNECTIONS
        [connections, list] = gen_phys_connections(striatum, attr);

        % TODO: FIX / UPDATE STATS GENERATION
        % gen_phys_connstats(striatum, connections)
    else 
        % GENERATE STATISTICAL INTRASTRIATAL CONNECTIONS
        [connections, list] = gen_stat_connections(attr);
    end
    
    % Generate striatal ID
    attr.id = strcat(datestr(now, 'yy.mm.dd_HH.MM'), '_', num2str(length(list.msn)), '+', num2str(length(list.fsi)));
    
    % Save striatum and connections 
    if attr.flags.save  
        attr.save_path = fullfile(attr.root, attr.id);
        mkdir(attr.save_path);

        if attr.flags.physical
            if attr.flags.progress
                fprintf('Saving striatum data… ')
                timer.save_str = tic;
            end

            % Save striatum data
            str_name = fullfile(attr.save_path, 'striatum.mat');
            save(str_name, 'striatum', '-v7.3');

            if attr.flags.progress
                fprintf('took %1.2f seconds.\n', toc(timer.save_str))
            end 
        end

        if attr.flags.progress
            fprintf('\nSaving connection data… ')
            timer.save_con = tic;
        end

        % Save neuron connection data
        con_name = fullfile(attr.save_path, 'connections.mat');
        save(con_name, 'connections', 'list', 'attr', '-v7.3');

        if attr.flags.progress
            fprintf('took %1.2f seconds. Done!\n', toc(timer.save_con))
        end
    end
    
else
    try
        % Load existing striatum and connection data
        fprintf('Loading existing data for striatum ID %s… ', s_ID);
        try
            load(fullfile(attr.root, s_ID, 'striatum.mat'));
        catch
            fprintf('no physical striatum found; trying to load statistical connections… ')
        end
             
        temp = attr.conn;
        
        load(fullfile(attr.root, s_ID, 'connections.mat'))
        fprintf('done!\n')
        
         attr.conn = temp;
        
    catch
        fprintf('failed!\n');
        error('Unable to load any existing striatal data');
    end
end

% GENERATE SPINECREATOR CONNECTION LISTS
if attr.flags.physical
    [connection_lists, attr] = gen_connection_lists(connections, list, attr, striatum); 
    if attr.flags.save
        % TODO: EXCLUDE THIS IF IT ALREADY EXISTS
        save_striatum(striatum, connection_lists, list, attr);
    end
else
    [connection_lists, ~] = gen_connection_lists(connections, list, attr);
end

if attr.flags.progress
    fprintf('Everything done! Total time elapsed: %1.2f minutes.\n', toc(timer.all) / 60)
end







% function attr = preserve(this)
% 
% end

% % Sanity checks
% if phys.size > 1000
%     error('Trying to create too large a striatum')
% elseif conn.bkg_msn > 100 || conn.bkg_fsi > 100
%     error('Trying to create too many background neurons')
% elseif phys.ch_width > 100 || phys.ch_width <= 0
%     error('Invalid channel width')
% else
%     
%     if ~exist('s_ID', 'var')
%         if flags.physical
%             % Create striatum structure
%             striatum = gen_phys_striatum(conn, phys, flags);
%             
%             if flags.progress
%                 timer.phys = tic;           % Timer for creation of connections
%                 fprintf('\nCreating physical connections…')
%                 fprintf('\nFor a full-size network (1mm³) this may take some time')
%                 fprintf('\nStart time: %s | ', datestr(now, 'HH:MM:SS'))
%             end
% 
%             %% GENERATE PHYSICAL INTRASTRIATAL CONNECTIONS
%             connections = gen_phys_connections(striatum, phys, conn, flags);
%             
%             % TODO: FIX / UPDATE STATS GENERATION
% %             gen_phys_connstats(striatum, connections)
%             
%             if flags.progress
%                 fprintf('\nCreating connections took %1.2f minutes', toc(timer.phys)/60 )
%             end
%         else
%             if flags.progress
%                 fprintf('\nCreating statistical connections…')
%                 timer.stat = tic;
%             end
%             
%             %% GENERATE STATISTICAL INTRASTRIATAL CONNECTIONS
%             connections = gen_stat_connections(stat, conn);
%             
%             if flags.progress
%                 fprintf('\nCreating connections took %1.2f seconds\n', toc(timer.stat))
%             end
%         end
%         
%         % Save connections to disk     
%         if flags.progress
%             timer.save = tic;
%             fprintf('\nSaving connection data… ')
%         end
%         
%         % FIXME: WHAT IS DIRNAME FOR STAT LIST
% %         filename = fullfile(attr.dirname, 'connections.mat');
% %         save(filename, 'connections');
%         
%         if flags.progress
%             fprintf('took %1.2f minutes. All done!\n', toc(timer.save)/60)
%         end
%     else
%         try
%             fprintf('Loading existing data for striatum ID %s… ', s_ID);
%             load(fullfile(attr.path, s_ID, 'striatum.mat'));
%             load(fullfile(attr.path, s_ID, 'connections.mat'))
%             fprintf('done!\n')
%         catch
%             fprintf('failed!\n');
%             error('Unable to load existing striatal data');
%         end
%         
% %         if flags.channel
% %             conn.ch_all = flags.channel;
% %             % TEST ONLY
% %             conn.bkg_msn = 50;
% %         end
%     end
% 
%     % TODO: Tidy up this iteration, make more flexible for future changes
%     
% %     % Create corticostriatal connections and output all connection lists
% %     if flags.density
% %         
% %         % Iterate active MSN density
% %         for m = 0:flags.density:(100 - flags.density)
% %             conn.bkg_msn = m;
% %             
% %             % Iterate active FSI density
% %             for f = 0:flags.density:100             
% %                 conn.bkg_fsi = f;
% %                 
% %                 if flags.width && conn.ch_all > 1
% %                     % Iterate channel width
% %                     for w = flags.width:flags.width:100
% %                         conn.ch_width = w;     
% %                         [connections, list] = gen_phys_connlist(striatum, connections, attr, conn, flags);                       
% %                     end
% %                 else
% %                     [connections, list] = gen_phys_connlist(striatum, connections, attr, conn, flags);
% %                 end
% %             end
% %         end
% %         
% %     elseif flags.width 
% %         % Iterate channel width
% %         for w = flags.width:flags.width:100          
% %             conn.ch_width = w;           
% %             [connections, list] = gen_phys_connlist(striatum, connections, attr, conn, flags);
% %         end
% %         
% %     else
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
%         %% WORKING HERE : SINGLE CONNECTION LIST
%         
%         list = get_neuron_list(connections);  
%         attr.id = strcat(datestr(now, 'yy.mm.dd_HH.MM'), '_', num2str(list.msn), '+', num2str(list.fsi));
%         
%         % Save striatum and connections 
%         if flags.save                 
%             if flags.progress
%                 fprintf('Saving striatum data… ')
%                 timer.save_str = tic;
%             end
% 
%             attr.save_path = fullfile(attr.root, attr.id);
%             mkdir(attr.save_path);
% 
%             filename = fullfile(attr.save_path, 'striatum.mat');
%             save(filename, 'striatum', 'phys');
% 
%             if flags.progress
%                 fprintf('took %1.2f seconds.\n', toc(timer.save_str))
%                 fprintf('\nSaving connection data… ')
%                 timer.save_con = tic;
%             end
% 
%             filename = fullfile(attr.save_path, 'connections.mat');
%             save(filename, 'connections');
% 
%             if flags.progress
%                 fprintf('took %1.2f minutes. Done!\n', toc(timer.save_con)/60)
%             end
%         end
%         
%         
%         connections = gen_connection_list(connections, list, conn);
%         
%         
% %         % Create single connection list with provided attributes
% %         [connections, list] = gen_phys_connlist(striatum, connections, attr, conn, flags);  
%         
%         
%         
%         
%         
%         
%         
%         
% %     end
%     
%     % Save all connections to disk
%     if flags.progress
%         fprintf('\nSaving complete connection details… ')
%     end
%     
%     if flags.save
%         save(fullfile(striatum.dirname, 'connections.mat'), 'connections', '-v7.3');
%     end
%     
%     if flags.progress
%         fprintf('done!\n')
%     end
%     
%     fprintf('All connection lists generated in %1.2f hours. Job complete.\n', toc(timer.all) / 3600)
% end