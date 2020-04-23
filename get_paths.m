function [s_dir, l_dir] = get_paths
% Generate paths to striatum and logfile data

HPC_NAME = 'ac1drb';
MAC_NAME = 'tacd';
LNX_NAME = 'dbuxton';

switch 7
    case exist(fullfile('/data', HPC_NAME), 'dir')
        % Iceberg / ShARC
        s_dir = fullfile('/data', HPC_NAME, 'striatums');

        switch 7
            case exist(fullfile('/fastdata-sharc', HPC_NAME), 'dir')
                % Iceberg
                l_dir = fullfile('/fastdata-sharc', HPC_NAME, 'output');

            case exist(fullfile('/fastdata', HPC_NAME), 'dir')
                % ShARC
                l_dir = fullfile('/fastdata', HPC_NAME, 'output');

            otherwise
                error('Unrecognised HPC system')
        end

    case exist('/Users', 'dir')
        % macOS
        s_dir = fullfile('/Users', MAC_NAME, 'Documents/Striatums');
        l_dir = fullfile('/Users', MAC_NAME, 'Documents/Logs/');

    case exist('/home', 'dir')
        % Linux desktop
        s_dir = fullfile('/home', LNX_NAME, 'Documents/Striatums');
        l_dir = fullfile('/home', LNX_NAME, 'Documents/Logs/');

    otherwise
        error('Unknown system or directory structure')
end