% Process discrete (spiking) log files produced by bulk experiments using
% physical striatal models

%% TODO LIST
% Add code to perform analysis on SGE
% Consider also using mtspecgrampt to get moving window of oscillations

%% PREAMBLE
% Reset initial state
clear variables; clc

% Add MATLAB path for HPC execution
addpath(genpath('~/MatLab'));
addpath(genpath('/home/ac1drb/MatLab'));

warning('off', 'MATLAB:MKDIR:DirectoryExists');

% Get path for loading striatum data and log files
[s_dir, l_dir] = get_paths;

%% CONFIGURATION
% Set basic parameters
attr.Striatum_ID = '20.04.10_17.00_84900+849';
attr.Experiment  = 'Physical';
attr.Channels    = 1;

% Striatum and log paths depend on current machine and experiment details
attr.Striatum_path = fullfile(s_dir, attr.Striatum_ID);
attr.Log_root = fullfile(l_dir, attr.Experiment);

% TODO: Move these to attr.
% Set key text for extraction of population name and trial information from logfile
csv_text = '_spike_';
trial_text = {'bkMSN', 'bkFSI'};
if attr.Channels > 1
    trial_text{end + 1} = 'wCH';
end

% Number of milliseconds over which to average spike counts
attr.BinWidth = 10;

% Set oscillation taper parameters
% W = 10;                       % Bandwidth
% T = 0.5;                      % Taper duration
% p = floor((2* W * T) - 1);    % Should be 2TW - 1
% params.tapers = [W, T, p];        
params.tapers = [5, 9];

% Duration (ms) of experiment to analyse for osicallations
params.time = 1000;

%% START
% Get all CSV files and remove non-log entries (e.g. connection lists)
csv_list = dir(fullfile(attr.Log_root, '/**/*.csv'));
csv_list(~contains({csv_list.folder}, '/log')) = [];

% Initialise output structures
results.Spikes = struct;
results.Oscillations = struct;

fprintf('Using striatum %s:\n', attr.Striatum_ID)

% Create logs structure if required
if ~isfield(attr, 'Logs')
    attr.Logs = struct;
end

% Get metadata for each log file
for i = 1:size(csv_list, 1)
    delims = strfind(csv_list(i).folder, '/');

    % Get experiment trial ID
    attr.Logs(i).Trial = strcat(...
        csv_list(i).folder(delims(end - 1) + 1 : delims(end) - 1));

    % Get path to each log file
    if ~strcmp(attr.Logs(i).Trial, attr.Experiment)        
        attr.Logs(i).Log_path = strcat(...
            csv_list(i).folder(delims(end - 1) + 1 : end), '/');
    else
        attr.Logs(i).Log_path = '/log/';
    end

    % Get log filename
    attr.Logs(i).Log_file = csv_list(i).name;  

    % Get population name
    attr.Logs(i).Population = attr.Logs(i).Log_file(...
        1 : strfind(attr.Logs(i).Log_file, csv_text) - 1);
    
%     for j = 1:length(trial_text)

    for j = 1:length(trial_text)
        % Get trial variables
        % From https://uk.mathworks.com/matlabcentral/answers/44049-extract-numbers-from-mixed-string

        % Get key text index
        idx = strfind(attr.Logs(i).Trial, trial_text{j});

        try
            % Extract number following key text
            attr.Logs(i).(trial_text{j}) = sscanf(...
                attr.Logs(i).Trial(idx(1) + length(trial_text{j}):end), '%g');
        catch
            % Ignore if single trial is being analysed
            fprintf('Unable to extract trial variable %s', trial_text{j});
        end
    end
end

% Process log files
for i = 1:size(attr.Logs, 2)
    
    % Load physical striatal data if different to previous trial
    if i == 1 || ~strcmp(attr.Logs(i).Trial, attr.Logs(i - 1).Trial)
        fprintf('Loading data for trial %s… ', attr.Logs(i).Trial)
        load(fullfile(attr.Striatum_path, 'neuron_data', attr.Logs(i).Trial, 'list.mat'), 'list');
        fprintf('done!\n')
    end
        
    % Reattribute and process spikes based on neuron ID in physical striatum  
    fprintf('    %s:\n', attr.Logs(i).Population)
    switch attr.Logs(i).Population
        case 'Striatum_D1'
            attr.Logs(i).Population_ID = 'd1';                    
        case 'Striatum_D2'
            attr.Logs(i).Population_ID = 'd2';
        case 'Striatum_FSI'
            attr.Logs(i).Population_ID = 'fsi';
        otherwise
            if contains(attr.Logs(i).Population, 'CH1_input')
                attr.Logs(i).Population = 'CH1_input';
            elseif contains(attr.Logs(i).Population, 'CH2_input')
                attr.Logs(i).Population = 'CH2_input';
            elseif contains(attr.Logs(i).Population, 'BKG')
                attr.Logs(i).Population = 'BKG_input';
            else
                error('Unknown striatal population encountered'); 
            end
    end
    
    % Load spike data and metadata
    [...
        raw_s, ...
        attr.Logs(i).Neurons_active, ...
        attr.Logs(i).Neurons_total, ...
        attr.Logs(i).Timestep, ...
        attr.Logs(i).Duration...
    ] ...
        = load_sc_discrete(fullfile(...
            attr.Log_root, ...
            attr.Logs(i).Log_path, ...
            attr.Logs(i).Log_file));
        
    % If the current logfile has no spikes, skip to the next logfile
    if isnan(attr.Logs(i).Neurons_total)
        fprintf('        (!) NO DATA FOUND (!)\n')
        continue
    end
            
    % If the current logfile is not an input
    if ~contains(attr.Logs(i).Population, 'input')
        % Reassign spike IDs based on physical striatum information
        [~, idx] = ismember(raw_s(:, 2), list.(attr.Logs(i).Population_ID)(:, 2));
        raw_s(:, 2) = list.(attr.Logs(i).Population_ID)(idx, 1);
    end
    
    % Create unique headers for this trial
    header = strcat(...
        attr.Experiment, '_', ...
        attr.Logs(i).Trial, '_', ...
        attr.Logs(i).Population);
    
    % If there are 2 or more channels and the current logfile is not an input
    if attr.Channels > 1 && ~contains(attr.Logs(i).Population, 'input')
        for c = 1:attr.Channels  
            ch = sprintf('ch%d', c);
            fprintf('        Processing channel %d: ', c)
            
            % Extract spikes from the current channel
            chn_s = raw_s(ismember(raw_s(:, 2), list.(ch).(attr.Logs(i).Population_ID)(:, 1)), :);
            
            % Create unique headers for this trial
%             header = strcat(header_prefix, '_', ch);
            [tot, avg, rol, spc, frq] = make_headers(strcat(header, '_', ch), attr.BinWidth);
            
            % Perform data analysis           
            fprintf('Spikes… ')
            [...
                results.Spikes.(tot), ...
                results.Spikes.(avg), ...
                results.Spikes.(rol)...
            ] ...
                = analyse_spikes(attr.Logs(i), chn_s, attr.BinWidth);
            fprintf('done! ')
            
            fprintf('Oscillations… ')
            [...
                results.Oscillations.(spc), ...
                results.Oscillations.(frq)...
            ]...
                = analyse_oscillations(attr.Logs(i), chn_s, params);
            fprintf('done!\n')
        end
    else
        fprintf('        Processing data: ')
        
        % Create unique headers for this trial
%         header = header_prefix;
        [tot, avg, rol, spc, frq] = make_headers(header, attr.BinWidth);
        
        fprintf('Spikes… ')
        [...
            results.Spikes.(tot), ...
            results.Spikes.(avg), ...
            results.Spikes.(rol)...
        ] ...
            = analyse_spikes(attr.Logs(i), raw_s, attr.BinWidth);
        fprintf('done! ')
        
        fprintf('Oscillations… ')
        [...
            results.Oscillations.(spc), ...
            results.Oscillations.(frq)...
        ]...
            = analyse_oscillations(attr.Logs(i), raw_s, params);
        fprintf('done!\n')
    end
    
    
    

    
    
    
    
      
    % TODO: Figure out why this isn't carving up spikes properly
%     for c = 1:attr.Channels       
% %         % If the current logfile is not an input
% %         if ~contains(attr.Logs(i).Population, 'input')
% %             ch = sprintf('ch%d', c);
% %             fprintf('        Processing channel %d: ', c)
% %             
% %             % Extract spikes from the current channel
% %             chn_s = raw_s(ismember(raw_s(:, 2), list.(ch).(attr.Logs(i).Population_ID)(:, 1)), :);
% %             
% %             % Create unique headers for this trial
% %             header = strcat(header, '_', ch);
% %         else 
% %             if c > 1
% %                 continue
% %             end
% %             fprintf('        Processing input: ')
% %             chn_s = raw_s; 
% %             
% % 
% %         end
%         
% %         % Spikes
% %         avg = strcat(header, '_SPmean');
% %         tot = strcat(header, '_SPtotal');       
% %         rol = strcat(header, sprintf('_SProll_%dms', attr.BinWidth));
% %         
% % %         % Oscillations
% % %         spc = strcat(header, '_OSspec');
% % %         frq = strcat(header, '_OSfreq');
% %         
% %         %% SPIKE ANALYSIS
% %         fprintf('Spikes… ')
% %                
% %         % Total spikes per millisecond
% %         results.Spikes.(tot) = [histcounts(raw_s(:,1), 0:attr.Logs(i).Duration)]'; %#ok<*NBRAK>
% % 
% %         % Mean spikes per neuron per second
% %         results.Spikes.(avg) = ...
% %             [histcounts(chn_s(:,1), 0:attr.Logs(i).Duration) ...
% %             / attr.Logs(i).Neurons_active * 1000]';
% % 
% %         % Rolling mean of spikes per neuron per 'binwidth' milliseconds
% %         results.Spikes.(rol) = movmean(results.Spikes.(avg), attr.BinWidth);
% % 
% %         fprintf('done! ')
%               
% %         %% OSCILLATION ANALYSIS
% %         fprintf('Oscillations… ')
% %         
% %         % Analyse only spikes from the last 'params.time' milliseconds
% %         trm_s = chn_s(chn_s(:,1) > max(chn_s(:,1)) - params.time, 1);
% % 
% %         [...
% %             results.Oscillations.(spc), ...
% %             results.Oscillations.(frq), ...
% %             ~...
% %         ] = mtspectrumpt(trm_s, params);
% % 
% %         % Calibrate spectral power to be accurate per *active* neuron
% %         results.Oscillations.(spc) = results.Oscillations.(spc) / attr.Logs(i).Neurons_active ^ 2;
% % 
% %         % Frequencies output must be transposed
% %         results.Oscillations.(frq) = results.Oscillations.(frq)';
% % 
% %         fprintf('done!\n')
%     end
end

% Append time to output structure
results.Spikes.(strcat(attr.Experiment, '_Time')) = (1:max([attr.Logs.Duration]))';

%% SAVE DATA
% Create directory for output data
out_dir = fullfile('~', attr.Experiment);
mkdir(out_dir);
fprintf('Saving data to CSV: ')

% Save spiking data
fprintf('Spikes… ')
struct2csv(results.Spikes, fullfile(...
    out_dir, strcat(attr.Experiment, '_spikes.csv')));
fprintf('done! ')

% Save oscillation data
fprintf('Oscillations… ')
struct2csv(results.Oscillations, fullfile(...
    out_dir, strcat(attr.Experiment, '_oscillations.csv')));
fprintf('done!\n')

function [tot, avg, rol, spc, frq] = make_headers(header, bw)
    % Spikes
    tot = strcat(header, '_SPtotal'); 
    avg = strcat(header, '_SPmean');          
    rol = strcat(header, sprintf('_SProll_%dms', bw));

    % Oscillations
    spc = strcat(header, '_OSspec');
    frq = strcat(header, '_OSfreq');
end

function [spk_tot, spk_avg, spk_rol] = analyse_spikes(log, spikes, bw)
    % Total spikes per millisecond
    spk_tot = [histcounts(spikes(:,1), 0:log.Duration)]'; %#ok<*NBRAK>

    % Mean spikes per neuron per second
    spk_avg = [histcounts(spikes(:,1), 0:log.Duration) / log.Neurons_active * 1000]';

    % Rolling mean of spikes per neuron per 'binwidth' milliseconds
    spk_rol = movmean(spk_avg, bw);
end

function [osc_spc, osc_frq] = analyse_oscillations(log, spikes, params)
     % Analyse only spikes from the last 'params.time' milliseconds
     % TODO: FIX THIS (analysing wrong part of ch1)
    trm_s = spikes(spikes(:,1) > max(spikes(:,1)) - params.time, 1);
    
    [osc_spc, osc_frq, ~] = mtspectrumpt(trm_s, params);
    
    % Calibrate spectral power to be accurate per *active* neuron
    osc_spc = osc_spc / log.Neurons_active ^ 2;

    % Frequencies output must be transposed
    osc_frq = osc_frq';
end