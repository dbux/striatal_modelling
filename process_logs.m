% Process discrete (i.e. spiking) log files produced by bulk experiments using
% physical striatal models

% % TODO: Add code to perform basic tasks only once on SGE
% 
% if isempty(getenv('SGE_TASK_ID'))
% %     error('Please run data analysis via Sun Grid Engine');
%     attr.SGE_ID = 95;   % TESTING ONLY
% else
%     attr.SGE_ID = str2double(getenv('SGE_TASK_ID')) - 1;
% end

% TODO: Consider also using mtspecgrampt to get moving window of
% oscillations

% Clear worksspace
clc; clear variables

% Add MATLAB path for ShARC execution
addpath(genpath('/home/ac1drb/MatLab'));

% Set basic parameters
attr.Striatum_ID = '20.02.25_10.51_84900+834_1CH_0sep_0overlap';
% attr.Striatum_path = strcat('/data/ac1drb/striatums/', attr.Striatum_ID, '/');
attr.Striatum_path = strcat('/Users/tacd/Documents/', attr.Striatum_ID, '/');
attr.Experiment = 'Physical_1CH';
attr.Channels = 1;

% Process log metadata
attr = get_metadata(attr);

% Process log files
[attr, results.Spikes, results.Oscillations] = process_data_spikes(attr);

% Create directory for output data
out_dir = strcat('~/', attr.Experiment);
mkdir(out_dir);
cd(out_dir);

% Save spiking data
fprintf('Saving spiking data to CSV? ')
output_name = strcat(attr.Experiment, '_spikes.csv');
struct2csv(results.Spikes, output_name);
fprintf('done!\n')

% Save oscillation data
fprintf('Saving oscillation data to CSV? ')
output_name = strcat(attr.Experiment, '_oscillations.csv');
struct2csv(results.Oscillations, output_name);
fprintf('done!\n')


% Extract extended attributes of bulk experiment log files
function attr = get_metadata(attr)

% Log path depends on current machine:
% ShARC
if exist('/fastdata/ac1drb', 'dir')
    attr.Log_root = strcat('/fastdata/ac1drb/output/', attr.Experiment, '/') ;
% Iceberg
elseif exist('/fastdata-sharc/ac1drb', 'dir')
    attr.Log_root = strcat('/fastdata-sharc/ac1drb/output/', attr.Experiment, '/');
% Local desktop (Linux)  
elseif exist('/home/dbuxton/', 'dir')
    attr.Striatum_path = strcat(...
        '/home/dbuxton/Dropbox/University/2019 (Sheffield) - Striatal oscillations/Striatums/', ...
        attr.Striatum_ID, '/');
    attr.Log_root = strcat('/home/dbuxton/Documents/Logs/', attr.Experiment, '/'); 
% Local desktop (Mac) 
elseif exist('/Users/tacd/', 'dir')
    attr.Striatum_path = strcat(...
        '/Users/tacd/Documents/Striatums/', ...
        attr.Striatum_ID, '/');
    attr.Log_root = strcat('/Users/tacd/Documents/Logs/', attr.Experiment, '/'); 
else
    error('Unable to determine logfile location');
end

% Get all CSV files and remove non-log entries (e.g. connection lists)
csv_list = dir(strcat(attr.Log_root, '/**/*.csv'));
csv_list(~contains({csv_list.folder}, '/log')) = [];

% Set key text for extraction of population name and trial information from logfile
csv_text = '_spike_';
trial_text = {'bkMSN', 'bkFSI'};

% Get logfile attributes
get_attributes(csv_list, csv_text);

    % Given a list of files and key text occuring immediately after the
    % population name, return each individual log file's attributes
    function get_attributes(file_list, pop_text)
        % Create logs structure if required
        if ~isfield(attr, 'Logs')
            attr.Logs = struct;
        end
        
        for i = 1:size(file_list, 1)
            delims = strfind(file_list(i).folder, '/');

            % Get experiment trial ID
            attr.Logs(i).Trial = strcat(...
                file_list(i).folder(delims(end - 1) + 1 : delims(end) - 1));

            % Get path to each log file
            if ~strcmp(attr.Logs(i).Trial , attr.Experiment)        
                attr.Logs(i).Log_path = strcat(...
                    file_list(i).folder(delims(end - 1) + 1 : end), '/');
            else
                attr.Logs(i).Log_path = '/log/';
            end

            % Get log filename
            attr.Logs(i).Log_file = file_list(i).name;  

            % Get population name
            attr.Logs(i).Population = attr.Logs(i).Log_file(...
                1 : strfind(attr.Logs(i).Log_file, pop_text) - 1); 
            
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
                    fprintf('Unable to extract trial variables, assuming single trial\n');
                end
            end
        end
    end
end


% Process spiking data
function[attr, output_spks, output_osc] = process_data_spikes(attr)
% TODO: Expand processing for each channel
% TODO: Get total and active neurons for each channel in multi-channel striatum
% TODO: Test spike processing with multi-channel data

% Load physical striatal data
fprintf('Loading striatal data from %s? ', attr.Striatum_path)
load(strcat(attr.Striatum_path, 'list.mat'), 'list');
fprintf('done!\n')

% Initialise output structures
output_spks = struct;
output_osc = struct;

% Process log files
for i = 1:size(attr.Logs, 2)
    % Get full path to log file
	logfile = char(strcat(...
        attr.Log_root, ...
        attr.Logs(i).Log_path, ...
        attr.Logs(i).Log_file));
    
    try
        % Load spiking data and associated metadata
        fprintf('Loading spiking data from %s? ', logfile)
        [spikes, ...
            attr.Logs(i).Neurons_active, ...
            attr.Logs(i).Neurons_total, ...
            attr.Logs(i).Timestep, ...
            attr.Logs(i).Duration] = load_sc_discrete(logfile);
        fprintf('done!\n')
        
        % Reattribute and process spikes based on neuron ID in physical striatum     
        switch attr.Logs(i).Population
            case 'Striatum_D1'
                attr.Logs(i).Population_ID = 'd1';                    
            case 'Striatum_D2'
                attr.Logs(i).Population_ID = 'd2';
            case 'Striatum_FSI'
                attr.Logs(i).Population_ID = 'fsi';
            otherwise
                error('Unknown population encountered');
        end
               
        process_spikes(attr.Logs(i).Population_ID);
        process_oscillations(spikes(:, 1));
    
    catch
        error('Unable to read data from %s', logfile)
    end
end

% Append time to output structure
output_spks.(strcat(attr.Experiment, '_', 'Time')) = (1:max([attr.Logs.Duration]))';


    % Reassign spike IDs and get detailed spiking data
    function data = process_spikes(msn)
        % Initialise data structure
        data = struct;
        
        % Reassign spike IDs based on physical striatum information
        [~, idx] = ismember(spikes(:, 2), list.(msn)(:, 2));
        spikes(:, 2) = list.(msn)(idx, 1);

        % Get mean and total spikes in each time bin
        if isfield(list, 'ch1')
            for c = 1:channels
                ch = sprintf('ch%d', c);
                
                [data.(ch).raw, data.(ch).avg, data.(ch).tot] = ...
                    mean_spikes(list.(ch).(msn)(:, 1));
            end
        else
            [data.raw, data.avg, data.tot] = mean_spikes(list.(msn)(:, 1));
        end
        
        % TODO: Move this into a separate nested funcion, update to allow
        % for multiple channels and different data
        % Create unique header text for this data
        header_spks = strcat(...
            attr.Experiment, '_', ...
            attr.Logs(i).Trial, '_', ...
            attr.Logs(i).Population, '_');
        
        % Append total and mean spiking data to output structure
        output_spks.(strcat(header_spks, 'spikes_total')) = data.tot;        
        output_spks.(strcat(header_spks, 'spikes_mean')) = data.avg;
    end


    % Extract raw, mean, and total number of spikes
    function [raw, avg, tot] = mean_spikes(n_list)
        % Number of milliseconds over which to average spike counts
        bw = 2;
    
        % Separate MSN spikes by channel into the spike struct
        raw = spikes(ismember(spikes(:, 2), n_list), :);
        
        % Preallocate avg array
        avg = zeros(attr.Logs(i).Duration, 1);
        tot = zeros(attr.Logs(i).Duration, 1);

        % Average number of spikes per second is calculated as:
        % 1. Total number of spikes in current time bin
        % 2. Divided by number of neurons in the current channel to get average spikes per neuron per time bin
        % 3. Multiplied by (1000 / bin width) to get average spikes per neuron per second       
        for t = 1:attr.Logs(i).Duration
            tot(t, 1) = length(find(...
                raw(:, 1) > (t - (bw / 2)) & ...
                raw(:, 1) < (t + (bw / 2))));
            
            avg(t, 1) = tot(t, 1) / attr.Logs(i).Neurons_active * (1000 / bw);
        end
    end

    function process_oscillations(spks)
        % Set taper parameters
%         W = 10;                       % Bandwidth
%         T = 0.5;                      % Taper duration
%         p = floor((2* W * T) - 1);    % Should be 2TW - 1
%         params.tapers = [W, T, p];        
        params.tapers = [5, 9];

        % Run multi-spectrum taper analysis on spike times
        [spec, freq, ~] = mtspectrumpt(spks, params);
        
        % Create unique header text for this data
        header_osc = strcat(...
            attr.Experiment, '_', ...
            attr.Logs(i).Trial, '_', ...
            attr.Logs(i).Population, '_');
        
        % Append total and mean spiking data to output structure
        output_osc.(strcat(header_osc, 'frequencies')) = freq';
        output_osc.(strcat(header_osc, 'spectrum')) = spec;
%         output_osc.(strcat(header_osc, 'rate')) = rate;    
    end
end