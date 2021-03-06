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

% Disable directory warnings
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% Get paths for striatum data and log files
[s_dir, l_dir] = get_paths;


%% CONFIGURATION
% Basic parameters
% attr.Striatum_ID = '20.04.10_17.00_84900+849';
% attr.Striatum_ID = '21.01.21_15.32_84900+811';

% attr.Experiment  = 'Physical_2CH_new';
% attr.Experiment  = '21.01.21_15.32_84900+811';

id.phys = '21.01.25_20.06_84900+845';
id.stat = '21.01.26_09.51_6000+60';

attr.Label       = 'PHYS';
attr.Channels    = 1;

if strfind(attr.Label, 'STAT')
    attr.Striatum_ID = id.stat;
    attr.Experiment  = id.stat;
elseif strfind(attr.Label, 'PHYS')
    attr.Striatum_ID = id.phys;
    attr.Experiment  = id.phys;
else
    error('Don''t know this label type')
end

% TODO: Extract number of channels from directory name

if attr.Channels == 1
    attr.C1_start = 1000;
    attr.C1_end   = 2000;
%     attr.C1_start = 0;
%     attr.C1_end   = 1000;
elseif attr.Channels == 2
    attr.C1_start = 500;
    attr.C1_end   = 2500;
    attr.C2_start = 1500;
    attr.C2_end   = 3500;
else
    fprintf('Unknown number of channels declared')
end

% TODO: Change trial_text and csv_text things into regular expressions

% Set key text for extraction of population name and trial information from logfile
attr.csv_text = '_spike_';
attr.trial_text = {'bkMSN', 'bkFSI'};
if attr.Channels > 1
    attr.trial_text{end + 1} = 'wCH';
end

% Analysis parameters
% Number of milliseconds over which to average spike counts
attr.BinWidth = 2;

% Oscillation taper parameters
% W = 10;                       % Bandwidth
% T = 0.5;                      % Taper duration
% p = floor((2 * W * T) - 1);	% Should be 2TW - 1
% params.tapers = [W, T, p];        
params.tapers = [5, 9];

% Striatum and log paths depend on current machine and experiment details
attr.Striatum_path = fullfile(s_dir, attr.Striatum_ID);
attr.Log_root = fullfile(l_dir, attr.Experiment);

% Duration (ms) of experiment to analyse for oscillations
% params.time = 1000;

%% GET TRIALS METADATA
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
        1 : strfind(attr.Logs(i).Log_file, attr.csv_text) - 1);

    for j = 1:length(attr.trial_text)
        % Get trial variables
        % From https://uk.mathworks.com/matlabcentral/answers/44049-extract-numbers-from-mixed-string

        % Get key text index
        idx = strfind(attr.Logs(i).Trial, attr.trial_text{j});

        try
            % Extract number following key text
            attr.Logs(i).(attr.trial_text{j}) = sscanf(...
                attr.Logs(i).Trial(idx(1) + length(attr.trial_text{j}):end), '%g');
        catch
            % Ignore if single trial is being analysed
            fprintf('Unable to extract trial variable %s', attr.trial_text{j});
        end
    end
end

% Process log files
for i = 1:size(attr.Logs, 2)
    
    % Load physical striatal data if different to previous trial
    if i == 1 || ~strcmp(attr.Logs(i).Trial, attr.Logs(i - 1).Trial)
        fprintf('Loading data for trial %s… ', attr.Logs(i).Trial)
%         load(fullfile(attr.Striatum_path, 'neuron_data', attr.Logs(i).Trial, 'list.mat'), 'list');
        load(fullfile(attr.Striatum_path, 'connections.mat'), 'list');
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
            % TEMP FOR extra input
%             elseif contains(attr.Logs(i).Population, 'MCtx_R2S')
%                 attr.Logs(i).Population = 'MCtx_input';
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
    
%     % Create unique headers for this trial
%     header = strcat(...
%         attr.Experiment, '_', ...
%         attr.Logs(i).Trial, '_', ...
%         attr.Logs(i).Population);
    % Create unique headers for this trial
%     header = strcat(attr.Striatum_ID, '_', attr.Experiment, '_', attr.Logs(i).Trial, '_', attr.Logs(i).Population);
    header = strcat(attr.Label, '_', attr.Logs(i).Trial, '_', attr.Logs(i).Population);
    
    % If there are 2 or more channels and the current logfile is not an input
    if attr.Channels > 1 && ~contains(attr.Logs(i).Population, 'input')
%     if attr.Channels > 1
        for c = 1:attr.Channels  
            ch = sprintf('ch%d', c);
            fprintf('        Processing channel %d: ', c)
            
            % Extract spikes from the current channel
            chn_s = raw_s(ismember(raw_s(:, 2), list.(ch).(attr.Logs(i).Population_ID)(:, 1)), :);
            
            % Create unique headers for this trial
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
            
%             plusone = results.Spikes.(avg)(results.Spikes.(avg) > 1)
% %             lastone = plusone(end) == plusone
%             idx = find(plusone(end) == plusone)
%             lastmainspike = idx(end)
            
            % TODO: Get last time index where avg spks/s is above 1, then make
            % that the final time index for CH1 oscillation analysis
            % OR come up with alternative method for analysing oscillation
            fprintf('Oscillations… ')
            if c == 1
                % TODO: Declare these up top somewhere
                solo_start = attr.C1_start;
                [solo_end, overlap_start] = deal(attr.C2_start);
                overlap_end = attr.C1_end;           
            elseif c == 2
                [solo_start, overlap_end] = deal(attr.C1_end);
                overlap_start = attr.C2_start;
                solo_end = attr.C2_end;
            else
                fprintf('Unknown channel for oscillation analysis')
            end
            
%             fprintf('Oscillations… ')
%             [...
%                 results.Oscillations.(spc), ...
%                 results.Oscillations.(frq)...
%             ]...
%                 = analyse_oscillations(attr.Logs(i), chn_s, params);
            fprintf('solo… ')
            
            % Append headers
            spc_solo = strcat(spc, '_solo');
            frq_solo = strcat(frq, '_solo');
            
            [...
                results.Oscillations.(spc_solo), ...
                results.Oscillations.(frq_solo)...
            ]...
            = analyse_oscillations(attr.Logs(i), chn_s, [solo_start, solo_end], params);

            fprintf('overlap… ')
            
            % Append headers
            spc_overlap = strcat(spc, '_overlap');
            frq_overlap = strcat(frq, '_overlap');
            
            [...
                results.Oscillations.(spc_overlap), ...
                results.Oscillations.(frq_overlap)...
            ]...
            = analyse_oscillations(attr.Logs(i), chn_s, [overlap_start, overlap_end], params);
        
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
        % Modify start and end time of oscillation analysis
        switch attr.Logs(i).Population
            case {'CH1_input', 'CH1A_input'}
%                 t_start = 1500;
%                 t_end   = 2500;
                t_start = 1500;
                t_end   = 1700;
            case 'CH2_input'
                t_start = 2500;
                t_end   = 3500;
            case 'BKG_input'
                t_start = 2500;
                t_end   = 3500;
            case 'MCtx_input'
                t_start = 1500;
                t_end   = 1700;
            otherwise
                t_start = attr.C1_start;
                t_end   = attr.C1_end;
        end
        [...
            results.Oscillations.(spc), ...
            results.Oscillations.(frq)...
        ]...
            = analyse_oscillations(attr.Logs(i), raw_s, [t_start, t_end], params);
        fprintf('done!\n')
    end
end

% Append time to output structure
% results.Spikes.(strcat(attr.Experiment, '_Time')) = (1:max([attr.Logs.Duration]))';
results.Spikes.Time = (1:max([attr.Logs.Duration]))';

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

%% FUNCTIONS

function [tot, avg, rol, spc, frq] = make_headers(header, bw)
    % Spikes
    tot = matlab.lang.makeValidName(strcat(header, '_SPtotal')); 
    avg = matlab.lang.makeValidName(strcat(header, '_SPmean'));          
    rol = matlab.lang.makeValidName(strcat(header, sprintf('_SProll_%dms', bw)));

    % Oscillations
    spc = matlab.lang.makeValidName(strcat(header, '_OSspec'));
    frq = matlab.lang.makeValidName(strcat(header, '_OSfreq'));
end

function [spk_tot, spk_avg, spk_rol] = analyse_spikes(log, spikes, bw)
    % Total spikes per millisecond
    spk_tot = [histcounts(spikes(:,1), 0:log.Duration)]'; %#ok<*NBRAK>

    % Mean spikes per neuron per second
    spk_avg = [histcounts(spikes(:,1), 0:log.Duration) / log.Neurons_active * 1000]';

    % Rolling mean of spikes per neuron per 'binwidth' milliseconds
    spk_rol = movmean(spk_avg, bw);
end

% function [osc_spc, osc_frq] = analyse_oscillations(log, spikes, params)
function [osc_spc, osc_frq] = analyse_oscillations(log, spikes, time, params)
     % Analyse only spikes from the last 'params.time' milliseconds   
%     trm_s = spikes(spikes(:,1) > max(spikes(:,1)) - params.time, 1);

    % Trim spikes to those in specified time region
    trm_s = spikes((time(2) > spikes(:,1) & spikes(:,1) > time(1)), 1);
    
    [osc_spc, osc_frq, ~] = mtspectrumpt(trm_s, params);
%     [osc_spc, osc_frq, ~] = mtspectrumpt(trm_s);
    
    % Calibrate spectral power to be accurate per *active* neuron
    osc_spc = osc_spc / log.Neurons_active ^ 2;

    % Frequencies output must be transposed
    osc_frq = osc_frq';
end