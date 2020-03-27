% Generates a model striatum and associated connections for either a
% physical or statistical version of the model

clear all; clc

addpath(genpath('~/MatLab'))

%clear striatum connections list
timer.all = tic;
rng('shuffle');

% Where to save striatum data files
if ~exist('~/Documents', 'dir')
    % HPC
    attr.path = '/data/ac1drb/striatums/';
else
    % Local
%     attr.path = '/Volumes/GoogleDrive/My Drive/1 - Projects/Striatal oscillations/Striatums/';
    attr.path = '~/Documents/Striatums/';
end

% Statistical attributes
% attr.num_msn = 6000;        % Number of MSNs
% attr.num_fsi = 60;          % Number of FSIs (should be 1% of MSNs)
attr.delay_min = 0.1;       % Minimum connection delay (must not be less than SpineCreator timestep)
attr.delay_mult = 0.4;      % Delay multiplier - should be less than 1

% Physical attributes
attr.size = 1000;            % Size of model striatum each side (?m) (Was 250?m in Humphries et al. 2009)
attr.min_dist = 10;         % Minimum distance between neurons (?m)
attr.centre_rad = 75;       % Radius of central region free from edge effects (?m)
attr.msn_density = 84900;   % Number of MSNs to place per mm^3 (should be 84,900)
attr.fsi_pct = 1;           % Percentage of MSNs to be added as FSIs  

% Connectivity attributes
attr.ch_seq = 1;            % Number of channels in the sequence
attr.ch_all = 1;            % Total number of channels (not including background)
attr.ch_inputs = 250;       % Number of input spike trains for each neuron
attr.ch_overlap = 0;        % Overlap (in % of striatal size) for two channels in physical model
attr.bkg_msn = 0;           % Percentage of MSNs to receive only background noise. Leave at 0 for no background.
attr.bkg_fsi = 0;           % Percentage of FSIs to receive only background noise. Leave at 0 for no background.
% attr.uni_bias = 100;        % Percentage bias towards SP connections being unidirectional
attr.max_bg = 500;          % Maximum number of connections allowed to a BG neuron (prevents max_spikes error in SpineCreator)
% attr.prn_src = 2;           % Prune SP connections going from SRC to DST channels
% attr.prn_dst = 1;
  
% Process flags
flags.phys_mod = 1;         % Create physical model?
flags.phys_ch = 0;          % Separate channels based on physical location? (TWO CHANNELS ONLY)
flags.debug = 0;            % Show detailed information during initialization?
flags.progress = 0;         % Show progress indicator? (Set to 0 for Iceberg)
flags.plot = 0;             % Show plot of striatum after generation?
flags.save = 1;             % Save connection lists to disk?
flags.binary = 1;           % Save binary versions of connection lists?
flags.density = 1;          % Create input lists for varying neural densities?


% Sanity checks
if flags.phys_ch && attr.ch_all ~= 2
    error('Trying to partition a physical model with more than two channels')
% elseif attr.ch_seq > attr.ch_all || attr.ch_seq < 2 || attr.ch_all < 2
%     error('Trying to create a model without enough channels')
elseif attr.size > 1000
    error('Trying to create too large a striatum')
elseif attr.bkg_msn > 100 || attr.bkg_fsi > 100
    error('Trying to create too many background neurons')
else
    % Proceed to creating a model
    if flags.phys_mod
        % % % % % % % % %
        % For PHYSICAL models
        % Create striatum structure
        striatum = gen_phys_striatum(attr, flags);

        % Create intrastriatal connections
        connections = gen_phys_connections(striatum, attr, flags);
        
        % Separate connections into channels and save connection lists
        if flags.density
            for i = 0:10:90
                for j = 0:10:90
                    % Iterate active neuron density
                    attr.bkg_msn = i;
                    attr.bkg_fsi = j;

%                     [connections, list] = gen_phys_connlist_1ch(striatum, connections, attr, flags);
                      [connections, list] = gen_phys_connlist(striatum, connections, attr, flags);
                end
            end
        else
            [connections, list] = gen_phys_connlist(striatum, connections, attr, flags);
        end
    else
        % % % % % % % % %
        % For STATISTICAL models
        % Create statistical connection lists
        connections = gen_stat_connections(attr, flags);
    end

    if flags.progress
        fprintf('Total time elapsed: %1.2f minutes.\n', toc(timer.all) / 60 )
    end
end


% % Optionally generate plot of E(c)
% % Page 9, figure 7
% figure(1); clf;
% semilogy(E(1,:), 'r', 'Linewidth', 2)
% xlabel('Distance d_s between somas (µm)')
% ylabel('Expected number of contacts')
% hold on
% grid on
% axis([0 600 0.0001 10])
% semilogy(E(2,:), 'b', 'Linewidth', 2)
% semilogy(E(3,:), 'g', 'Linewidth', 2)
% semilogy(E(4,:), 'k', 'Linewidth', 2)
% legend('MSN-MSN', 'FSI-MSN', 'FSI-FSI', 'FSI gap')


% Plot the striatum if desired
if flags.plot == 1
    figure(2); clf;
    scatter3(striatum.neurons(:,1),striatum.neurons(:,2),striatum.neurons(:,3), 0.1)
    hold on
    scatter3(striatum.neurons_centre(:,1),striatum.neurons_centre(:,2),striatum.neurons_centre(:,3), 20, 'o', 'r', ...
      'MarkerFaceColor', 'r')
    scatter3(striatum.neurons(striatum.linear_centre(striatum.linear_centre(:,2)==fsi),1),...
        striatum.neurons(striatum.linear_centre(striatum.linear_centre(:,2)==fsi),2),...
        striatum.neurons(striatum.linear_centre(striatum.linear_centre(:,2)==fsi),3), 20, 'o', 'm', 'MarkerFaceColor', 'm')

    res = 25;
    r = 75 * ones(res, res);
    [th, phi] = meshgrid(linspace(0, 2*pi, res), linspace(-pi, pi, res));
    [x,y,z] = sph2cart(th, phi, r);
    x = x + strinf.centre;
    y = y + strinf.centre;
    z = z + strinf.centre;
    surface(x,y,z,'FaceColor', 'none')   
end