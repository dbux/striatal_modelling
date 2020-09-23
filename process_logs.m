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
attr.Striatum_ID = '20.04.10_17.00_84900+849_2CH';
attr.Experiment  = 'Physical_2CH';
attr.Channels    = 2;

% Striatum and log paths depend on current machine and experiment details
attr.Striatum_path = fullfile(s_dir, attr.Striatum_ID);
attr.Log_root = fullfile(l_dir, attr.Experiment);

% TODO: Move these to attr.
% Set key text for extraction of population name and trial information from logfile
csv_text = '_spike_';
trial_text = {'bkMSN', 'bkFSI', 'wCH'};

% Number of milliseconds over which to average spike counts
bw = 10;

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
            error('Unknown striatal population encountered');
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
        
    % Reassign spike IDs based on physical striatum information
    [~, idx] = ismember(raw_s(:, 2), list.(attr.Logs(i).Population_ID)(:, 2));
    raw_s(:, 2) = list.(attr.Logs(i).Population_ID)(idx, 1);
        
    for c = 1:attr.Channels
        ch = sprintf('ch%d', c);
        fprintf('        Processing channel %d: ', c)

        % Create unique headers for this trial
        header = strcat(...
            attr.Experiment, '_', ...
            attr.Logs(i).Trial, '_', ...
            attr.Logs(i).Population, '_', ...
            ch);

        % Spikes
        avg = strcat(header, '_SPmean');
        tot = strcat(header, '_SPtotal');       
        rol = strcat(header, sprintf('_SProll_%dms', bw));
        
        % Oscillations
        spc = strcat(header, '_OSspec');
        frq = strcat(header, '_OSfreq');
        
        %% SPIKE ANALYSIS
        fprintf('Spikes… ')

        % Extract spikes from the current channel
        chn_s = raw_s(ismember(raw_s(:, 2), list.(ch).(attr.Logs(i).Population_ID)(:, 1)), :);

        % Total spikes per millisecond
        results.Spikes.(tot) = [histcounts(chn_s(:,1), 0:attr.Logs(i).Duration)]'; %#ok<*NBRAK>

        % Mean spikes per neuron per second
        results.Spikes.(avg) = ...
            [histcounts(chn_s(:,1), 0:attr.Logs(i).Duration) ...
            / attr.Logs(i).Neurons_active * 1000]';
        
        % Rolling mean of spikes per neuron per 'binwidth' milliseconds
        results.Spikes.(rol) = movmean(results.Spikes.(avg), bw);
        
        fprintf('done! ')
              
        %% OSCILLATION ANALYSIS
        fprintf('Oscillations… ')
        % Analyse only spikes from the last 'params.time' milliseconds
        trm_s = chn_s(chn_s(:,1) > max(chn_s(:,1)) - params.time, 1);

        [...
            results.Oscillations.(spc), ...
            results.Oscillations.(frq), ...
            ~...
        ] = mtspectrumpt(trm_s, params);
    
        % Calibrate spectral power to be accurate per *active* neuron
        results.Oscillations.(spc) = results.Oscillations.(spc) / attr.Logs(i).Neurons_active ^ 2;
                     
        % Frequencies output must be transposed
        results.Oscillations.(frq) = results.Oscillations.(frq)';
        
        fprintf('done!\n')
    end
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