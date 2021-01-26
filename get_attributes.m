%% Striatal generation configuration
% General attributes
[attr.root, ~] = get_paths;     % Get path for saving striatum data files

%% Physical striatum attributes
attr.phys.size        = 1000;   % Size of model striatum each side (μm) (Was 250μm in Humphries et al. 2009)
attr.phys.min_dist    = 10;   	% Minimum distance between neurons (μm)
attr.phys.centre_rad  = 75;     % Radius of central region free from edge effects (μm)
attr.phys.msn_density = 84900;  % Number of MSNs to place per mm^3 (should be 84,900)
attr.phys.fsi_pct     = 1;      % Percentage of MSNs to be added as FSIs  
attr.phys.ch_width    = 50;     % Percentage width of each channel in two-channel physical model (0 means no channel inputs).
attr.phys.cv_msnmsn   = 1;      % MSN-MSN conductance velocity (m/s)
attr.phys.cv_fsimsn   = 1;      % FSI-MSN conductance velocity (m/s)    
attr.phys.cv_fsifsi   = 1;      % FSI-FSI conductance velocity (m/s)

%% Statistical striatum attributes
attr.stat.num_msn = 6000;       % Number of MSNs
attr.stat.num_fsi = 60;         % Number of FSIs (should be 1% of MSNs)

attr.stat.con_msnmsn = 728;     % Expected number of connections of each type  
attr.stat.con_fsimsn = 3017;    % Values from Humphries, Wood and Gurney (2010), page 11, table 5
attr.stat.con_fsifsi = 12.8;
attr.stat.con_fsigap = 0.65;

%% Connectivity attributes
attr.conn.ch_all      = 1;      % Total number of input channels (not including background)
attr.conn.bkg_msn     = 40;     % Percentage of MSNs to receive only background noise. Leave at 0 for no background.
attr.conn.bkg_fsi     = 0;      % Percentage of FSIs to receive only background noise. Leave at 0 for no background.
attr.conn.delay_min   = 0.1;    % Minimum connection delay (must not be less than SpineCreator timestep)
attr.conn.delay_mult  = 1;      % Connection list delay factor

%% Process flags
attr.flags.physical   = 1;      % Create physical-based striatum?
attr.flags.phys_ch    = 1;      % Physically partition two-channel striatum?
attr.flags.debug      = 0;      % Show detailed information during initialization?
attr.flags.progress   = 1;      % Show progress indicator? (Set to 0 for Iceberg)
attr.flags.save       = 1;      % Save connection lists to disk?
attr.flags.binary     = 1;      % Save binary versions of connection lists?
attr.flags.density    = 0;      % Vary neural density by this interval (0 for single list)
attr.flags.width      = 0;      % Vary channel with by this interval (0 for single list)