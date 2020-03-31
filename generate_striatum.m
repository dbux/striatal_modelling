% Generates a 'physical' striatal microcircuit and associated connection 
% lists based on Humphries, Wood & Gurney (2010)

% Reset initial state
clear variables; clc
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
 
% % % % % % % % % %
% DEFINE ATTRIBUTES

% Physical striatum attributes
attr.size        = 500;    % Size of model striatum each side (μm) (Was 250μm in Humphries et al. 2009)
attr.min_dist    = 10;      % Minimum distance between neurons (μm)
attr.centre_rad  = 75;      % Radius of central region free from edge effects (μm)
attr.msn_density = 84900;   % Number of MSNs to place per mm^3 (should be 84,900)
attr.fsi_pct     = 1;       % Percentage of MSNs to be added as FSIs  

% Connectivity attributes
% attr.ch_seq     = 1;          % Number of channels in the sequence
attr.ch_all     = 1;        % Total number of channels (not including background)
% attr.ch_inputs  = 250;    	% Number of input spike trains for each neuron
attr.ch_overlap = 0;        % Overlap (in % of striatal size) for two channels in physical model
attr.bkg_msn    = 0;        % Percentage of MSNs to receive only background noise. Leave at 0 for no background.
attr.bkg_fsi    = 0;        % Percentage of FSIs to receive only background noise. Leave at 0 for no background.
% attr.uni_bias  = 100;         % Percentage bias towards SP connections being unidirectional
attr.max_bg     = 500;      % Maximum number of connections allowed to a BG neuron (prevents max_spikes error in SpineCreator)

% Process flags
flags.phys_ch   = 0;        % Separate channels based on physical location? (TWO CHANNELS ONLY)
flags.debug     = 0;        % Show detailed information during initialization?
flags.progress  = 1;        % Show progress indicator? (Set to 0 for Iceberg)
flags.plot      = 0;        % Show plot of striatum after generation?
flags.save      = 1;        % Save connection lists to disk?
flags.binary    = 1;        % Save binary versions of connection lists?
flags.density   = 1;        % Create input lists for varying neural densities?

% % % % % % % % % %
% BEGIN

% Sanity checks
if flags.phys_ch && attr.ch_all ~= 2
    error('Trying to partition a physical model without two channels')
elseif attr.size > 1000
    error('Trying to create too large a striatum')
elseif attr.bkg_msn > 100 || attr.bkg_fsi > 100
    error('Trying to create too many background neurons')
else
    % Create striatum structure
    striatum = gen_phys_striatum(attr, flags);
    
    % Create intrastriatal connections
    connections = gen_phys_connections(striatum, attr, flags);
    
end

% % % % % % % % % %
% FUNCTIONS

function[striatum] = gen_phys_striatum(attr, flags)
    % Generates a model striatum with topology based on the description in
    % page 7 of Humphries, Wood & Gurney (2010)
    %
    % !IMPORTANT! This program requires the 'Geom3D' toolbox available from:
    % http://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
    % Specifically, the files 'distancePoints3d.m' and 'vectorNorm3d.m' must be
    % available for calculating neuron distances.

    timer.str = tic;                                                % Timer for striatum creation
    reverseStr = '';                                                % Required string for progress indicator

    pop.msn = floor(round(attr.size^3 * attr.msn_density) / 1e9);   % Number of MSNs to place
    while mod((pop.msn/2), attr.ch_all) ~= 0                        % Ensure the number of MSNs will fit nicely into the number of channels
        pop.msn = pop.msn + 1;
    end

    pop.max = round(pop.msn + (pop.msn/10));                        % Maximum possible population (for matrix preallocation only)
    strinf.centre = attr.size/2;                                    % Centre point of striatum
    strinf.id = datestr(now, 'yy.mm.dd_HH.MM');                     % Unique ID for this striatum

    % How neurons will be represented numerically
    msn = 1;                
    fsi = 3;

    % Preallocate striatal array - initialize as uint8 to save memory and time
    % Striatal centre arrays small enough to not require preallocation
    striatum.main = zeros(attr.size,attr.size,attr.size, 'uint16');
    striatum.linear = zeros(pop.max,1, 'uint8');
    striatum.linear_centre = [];
    striatum.neurons = zeros(pop.max,3, 'uint16');
    striatum.neurons_centre = [];
    strinf.ptr = 1;                 % Points to the next available striatum entry

    if flags.progress
        fprintf('\nInitializing %dμm³ striatum with %d MSNs + %d%% FSIs… ', attr.size, pop.msn, attr.fsi_pct)
    end

    % Place each neuron individually
    for i = 1:pop.msn  
        [striatum, strinf] = gen_place_neuron(striatum, strinf, attr, msn);

        % Roll the dice to see if we will also be generating an FSI; FSI
        % creation procedure is the same as above
        if randi(100) <= attr.fsi_pct
            [striatum, strinf] = gen_place_neuron(striatum, strinf, attr, fsi);       
        end

        % Display percentage complete
        % http://undocumentedmatlab.com/blog/command-window-text-manipulation
        if flags.debug == 0 && flags.progress      
            percentDone = 100 * i / pop.msn;
            msg = sprintf('percent done: %3.1f', percentDone);
            fprintf([reverseStr, msg])
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end

    % Remove any empty space at the end of the neuron list and location register
    striatum.linear = striatum.linear(striatum.linear~=0);
    striatum.neurons = striatum.neurons(any(striatum.neurons,2),:);

    % Get some additional neuron population info
    pop.fsi = length(find(striatum.linear==fsi));
    pop.fsi_centre = length(find(striatum.linear_centre(:,2)==fsi));
    pop.all = length(striatum.linear);
    pop.all_centre = length(striatum.linear_centre);

    if flags.progress
        fprintf('\nStriatum initialization took %1.2f minutes - %d of %d neurons are FSIs', toc(timer.str)/60, pop.fsi, pop.all)
        fprintf('\nA total of %d neurons are in the centre region, of which %d are FSIs\n', pop.all_centre, pop.fsi_centre)
    end

    % Save striatum to disk
    timer.save = tic;
    if flags.progress
        fprintf('Saving striatum data… ')
    end
    striatum.dirname = [attr.path num2str(strinf.id) '_' ...
        num2str(pop.msn) '+' num2str(pop.fsi) '_' num2str(attr.ch_all) 'CH_' num2str(flags.phys_ch) 'sep_' num2str(attr.ch_overlap) 'overlap/'];
    mkdir(striatum.dirname);
    filename = [striatum.dirname '/striatum.mat'];
    save(filename, 'striatum', 'strinf', 'attr');
    if flags.progress
        fprintf('took %1.2f seconds. Done!\n', toc(timer.save))
    end
end

% Only called by gen_phys_striatum but runs slowly as a nested function    
function[striatum, strinf] = gen_place_neuron(striatum, strinf, attr, n_type)
    % Given a striatum structure and associated information, will place a
    % neuron in the physical striatum ensuring a minimum distance from all
    % other neurons and return the co-ordinates of that neuron

    % Reset 'neuron OK' flag before entering while-loop
    n_ok = 0;

    while n_ok == 0
        % Generate random xyz co-ordinates for neuron placement
        n_loc = randi(attr.size,1,3);
        % Assume that the co-ords are OK unless we decide otherwise
        n_ok = 1;

        % If those co-ords are already occupied, generate a new set
        while striatum.main(n_loc(1),n_loc(2),n_loc(3)) ~= 0
            n_loc = randi(attr.size,1,3);
        end

        % Create a set of 'buffer' co-ords that extend the minimum
        % separation distance away from the potential neuron co-ords in all directions
        for j = 1:3
            if n_loc(j) - attr.min_dist < 1
                buffer{j}(1) = 1;
            else
                buffer{j}(1) = n_loc(j) - attr.min_dist;
            end

            if n_loc(j) + attr.min_dist > attr.size
                buffer{j}(2) = attr.size;
            else
                buffer{j}(2) = n_loc(j) + attr.min_dist;
            end
        end

        % Test to see if any neuron already in this area would be too close
        % to the new neuron. Test is for distance, not just presence,
        % because diagonal distances are greater and may be acceptable
        for x = buffer{1}(1):buffer{1}(2)
            for y = buffer{2}(1):buffer{2}(2)
                for z = buffer{3}(1):buffer{3}(2)
                    if striatum.main(x,y,z) ~= 0
                        if distancePoints3d([n_loc(1) n_loc(2) n_loc(3)], [x y z]) < attr.min_dist
    %                         if flags.debug == 1
    %                             fprintf('\nNeuron %d at [%d, %d, %d] only %1.2fµm from neuron at [%d, %d, %d] - trying again', ...
    %                                 i, n(1), n(2), n(3), distancePoints3d([n(1) n(2) n(3)], [x y z]), x, y, z)
    %                         end

                            % If minimum distance is violated then set
                            % 'neuron not OK' and start again
                            n_ok = 0;
                        end
                    end
                end
            end
        end
    end

    % If the new neuron is close to the centre of the striatum, add it to
    % the list of central neurons
    if distancePoints3d([strinf.centre strinf.centre strinf.centre], [n_loc(1) n_loc(2) n_loc(3)]) <= attr.centre_rad
        striatum.linear_centre = [striatum.linear_centre; strinf.ptr, n_type];
        striatum.neurons_centre = [striatum.neurons_centre ; n_loc(1) n_loc(2) n_loc(3)];
    end

    % Everything checks out! Add the neuron to the striatum and the neuron
    % location register
    striatum.main(n_loc(1),n_loc(2),n_loc(3)) = n_type;
    striatum.linear(strinf.ptr) = n_type;
    striatum.neurons(strinf.ptr,:) = [n_loc(1) n_loc(2) n_loc(3)];
    strinf.ptr = strinf.ptr + 1;
end

function[connections] = gen_phys_connections(striatum, attr, flags)
    %  Humphries, Wood & Gurney (2010)
    %  Page 10

    warning('off', 'stats:pdist2:DataConversion');

    timer.conn = tic;           % Timer for creation of connections
    reverseStr = '';            % Required string for progress indicator
    E = gen_e;                  % Generate values for E
    ptr = 1;                    % Points to the next available synaptic connection entry                              
    msn = 1;                    % How neurons will be represented numerically            
    fsi = 3;

    delay_mult = 0.02;      % Delay multiplier - multiply connection distances by this for delay value
    delay_min = 0.1;        % Minimum delay (must be at least as large as SpineCreator timestep value)
    r = rand;               % Use one random number for all delays so they are in relative proportion

    % Prellocate synaptic connection lists
    connections.msnmsn = zeros((1000*length(striatum.linear)),3);
    connections.fsimsn = zeros((1000*length(striatum.linear)),3);
    connections.fsifsi = zeros((1000*length(striatum.linear)),3);
    connections.gap = [];       % Gap connection list does not need preallocation

    if flags.progress
        fprintf('\nCreating connections…')
        fprintf('\n(!) For a full-size network (1mm³) creating connections may take 1-3 hours (!)')
        fprintf('\nStart time: %s | ', datestr(now, 'HH:MM:SS'))
    %     fprintf('\nEstimated completion time: %1.2f - %1.2f minutes... ', (0.09*length(striatum.linear))/60, (0.15*length(striatum.linear))/60)
    end

    % For each neuron...
    for i = 1:length(striatum.linear)
        timer.neuron = tic;

        coninf.count = 0;       % Counts number of connections created 
        coninf.dist_all = [];   % Distance from this neuron to all other neurons
        coninf.exp_all = [];    % Expected number of connections of all types to all other neurons 
        coninf.exp_this = [];   % Expected number of connections of appropriate types to all other neurons
        coninf.trials = [];     % Number of connection trials to perform for all potential connections
        coninf.prob = [];       % Probability of making a connection at each trial to all other neurons

        if flags.debug
            fprintf('\nConnecting neuron %d of %d… ', i, length(striatum.linear))
        end

        % Get the distance from current neuron to all other neurons. Must be a
        % round number to look up expected connections from E later.    
    %     coninf.dist_all = round(distancePoints3d(striatum.neurons(i,:), striatum.neurons));
        coninf.dist_all = pdist2(striatum.neurons(i,:), striatum.neurons);       

        % This is a bit of a hack - to prevent errors from unwanted 0-distance
        % connections (i.e. self-connections) we set the distance to 2000
        % (maximum distance, with 0 expected connections for all types)
        coninf.dist_all(coninf.dist_all==0) = length(E);

        % Get the type of connection to all other neurons
        % MSN = 1, FSI = 3
        % So MSN-MSN = 1+1 = 2; FSI-MSN = 3+1 = 4; FSI-FSI = 3+3 = 6
        % (Must be converted to UINT32 to avoid forcing connection lists into UINT8 format)
        coninf.type = uint32(striatum.linear(i) + striatum.linear);

        % If current neuron is 1 (MSN), remove connections of type 4 (MSN-FSI) as they are invalid
        if striatum.linear(i) == msn
            coninf.type(coninf.type==(msn+fsi)) = 0;
        end 

        % Get the expected number of connections from this neuron to every
        % other neuron based on their relative distance. At this stage we don't
        % discriminate based on the type of connection being made, and obtain
        % expected connection values for all connection types.
        coninf.exp_all = E(:,round(coninf.dist_all))';

        % Get the values of only the desired connection types. For FSI-FSI 
        % connections we want values for both synapses and gap junctions so two 
        % columns are needed
        coninf.exp_this(coninf.type==(msn+msn),1) = coninf.exp_all(coninf.type==(msn+msn),1);   
        coninf.exp_this(coninf.type==(fsi+msn),1) = coninf.exp_all(coninf.type==(fsi+msn),2);   
        coninf.exp_this(coninf.type==(fsi+fsi),1) = coninf.exp_all(coninf.type==(fsi+fsi),3);   
        coninf.exp_this(coninf.type==(fsi+fsi),2) = coninf.exp_all(coninf.type==(fsi+fsi),4);

        % Because some expected connection values are >1, a standard method for
        % generating probabilistic connections is required:
        % Let n be the expected number of connections between two neurons.
        % (2*ceiling-n) trials are run, each with a (0.5 * n/ceiling-n) chance 
        % of making a connection.
        coninf.trials = 2 * ceil(coninf.exp_this);
        coninf.prob = 0.5 * (coninf.exp_this ./ ceil(coninf.exp_this));

        % Run the trials for each pair and generate connections. Store the
        % connection pair, type, and the distance between the two neurons
        for j = 1:max(coninf.trials(:,1))
            for k = 1:length(coninf.trials)
                if coninf.trials(k,1) ~= 0
                    if rand < coninf.prob(k,1)
                        coninf.dist_this = coninf.dist_all(k);
                        if coninf.type(k) == (msn+msn)
                            connections.msnmsn(ptr,:) = [i k coninf.dist_this];
                        elseif coninf.type(k) == (fsi+msn)
                            connections.fsimsn(ptr,:) = [i k coninf.dist_this];
                        elseif coninf.type(k) == (fsi+fsi)
                            connections.fsifsi(ptr,:) = [i k coninf.dist_this];
                        else
                            error('Unexpected connection type!')
                        end
                        ptr = ptr+1;
                        coninf.count = coninf.count + 1;
                    end
                    coninf.trials(k,1) = coninf.trials(k,1) - 1;
                end
            end
        end

        % If the current neuron is an FSI, run connection trials for gap junctions
        if striatum.linear(i) == fsi
            for j = 1:max(coninf.trials(:,2))
                for k = 1:length(coninf.trials)
                    if coninf.trials(k,2) ~= 0
                        if rand < coninf.prob(k,2)
                            coninf.dist_this = coninf.dist_all(k);
                            connections.gap = [connections.gap; i k coninf.dist_this];
                            coninf.count = coninf.count + 1;
                        end
                        coninf.trials(k,2) = coninf.trials(k,2) - 1;
                    end
                end
            end
        end

        if flags.debug
            fprintf('made %d connections in %1.2f seconds', coninf.count, toc(timer.neuron))
        end

        if flags.debug == 0 && flags.progress
            percentDone = 100 * i / length(striatum.linear);
            msg = sprintf('Percent done: %3.1f', percentDone);
            fprintf([reverseStr, msg])
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end

    % Remove empty space at the end of connection lists
    connections.msnmsn = connections.msnmsn(any(connections.msnmsn,2),:);
    connections.fsimsn = connections.fsimsn(any(connections.fsimsn,2),:);
    connections.fsifsi = connections.fsifsi(any(connections.fsifsi,2),:);

    % Convert connections distances to delay
    connections.msnmsn(:,3) = connections.msnmsn(:,3) .* delay_mult .* r + delay_min;
    connections.fsimsn(:,3) = connections.fsimsn(:,3) .* delay_mult .* r + delay_min;
    connections.fsifsi(:,3) = connections.fsifsi(:,3) .* delay_mult .* r + delay_min;

    if flags.progress
        fprintf('\nCreating connections took %1.2f minutes', toc(timer.conn)/60 )
    end

    % Save connections to disk
    timer.save = tic;
    if flags.progress
        fprintf('\nSaving connection data… ')
    end
    filename = [striatum.dirname '/connections.mat'];
    save(filename, 'connections');
    if flags.progress
        fprintf('took %1.2f minutes. All done!\n', toc(timer.save)/60)
    end

    % % % % % % % % %
    % Generate connectivity statistics
    statname = [striatum.dirname '/connection_stats.csv'];
    fc = fopen(statname, 'at+');
    fprintf(fc, '\n%d MSNs, %d FSIs, FSI ratio %d%%', ...
        length(find(striatum.linear==msn)), length(find(striatum.linear==fsi)), attr.fsi_pct);
    fprintf(fc, '\nConnection type, No. contacts, STD, Distance (Âµm), STD');

    % MSNs - 1 MSN
    try
        msn1msn = gen_conn_stats(connections.msnmsn, striatum, 0);   
        fprintf(fc, '\nMSNs - 1 MSN, %1.2f, %1.2f, %1.2f, %1.2f', mean(msn1msn.numbers), std(msn1msn.numbers), ...
            mean(msn1msn.dists), std(msn1msn.dists));
    catch
        fprintf(fc, '\nError generating MSNs - 1 MSN connection statistics!');
    end

    % FSIs - 1 MSN
    try
        fsi1msn = gen_conn_stats(connections.fsimsn, striatum, 0);
        fprintf(fc, '\nFSIs - 1 MSN, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsi1msn.numbers), std(fsi1msn.numbers), ...
            mean(fsi1msn.dists), std(fsi1msn.dists));
    catch
        fprintf(fc, '\nError generating FSIs - 1 MSN connection statistics!');
    end

    % 1 FSI - MSNs
    try
        fsimsns = gen_conn_stats(connections.fsimsn, striatum, 1);
        fprintf(fc, '\n1 FSI - MSNs, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsimsns.numbers), std(fsimsns.numbers), ...
            mean(fsimsns.dists), std(fsimsns.dists));
    catch
        fprintf(fc, '\nError generating 1 FSI - MSNs connection statistics!');
    end

    % FSIs - 1 FSI
    try
        fsi1fsi = gen_conn_stats(connections.fsifsi, striatum, 0);
        fprintf(fc, '\nFSIs - 1 FSI, %1.2f, %1.2f, %1.2f, %1.2f', mean(fsi1fsi.numbers), std(fsi1fsi.numbers), ...
            mean(fsi1fsi.dists), std(fsi1fsi.dists));
    catch
        fprintf(fc, '\nError generating FSIs - 1 FSI connection statistics!');
    end

    % FSI gap
    try
        fsigap = gen_conn_stats(connections.gap, striatum, 0);
        fprintf(fc, '\nFSI gap, %1.2f, %1.2f, %1.2f, %1.2f\n', mean(fsigap.numbers), std(fsigap.numbers), ...
            mean(fsigap.dists), std(fsigap.dists));
    catch
        fprintf(fc, '\nError generating FSI gap connection statistics!');
    end

    fclose(fc);

    % DISTANCES ARE WRONG

    if flags.progress
        fprintf('\nStriatal connectivity statistics (%d%% FSI):', attr.fsi_pct)
        fprintf('\n----------------------------------------------------------')
        fprintf('\nType\t\t| No. of contacts\t| Distance (µm)')
        fprintf('\n----------------------------------------------------------\n')
        fprintf('MSNs - 1 MSN\t| %1.2f, ± %1.2f\t| %1.2f ± %1.2f\n', mean(msn1msn.numbers), std(msn1msn.numbers), ...
            mean(msn1msn.dists), std(msn1msn.dists))
        fprintf('FSIs - 1 MSN\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsi1msn.numbers), std(fsi1msn.numbers), ...
            mean(fsi1msn.dists), std(fsi1msn.dists))
        fprintf('1 FSI - MSNs\t| %1.2f, ± %1.2f\t| %1.2f ± %1.2f\n', mean(fsimsns.numbers), std(fsimsns.numbers), ...
            mean(fsimsns.dists), std(fsimsns.dists))
        fprintf('FSIs - 1 FSI\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsi1fsi.numbers), std(fsi1fsi.numbers), ...
            mean(fsi1fsi.dists), std(fsi1fsi.dists))
        try
            fprintf('FSI gap junc\t| %1.2f, ± %1.2f\t\t| %1.2f ± %1.2f\n', mean(fsigap.numbers), std(fsigap.numbers), ...
            mean(fsigap.dists), std(fsigap.dists))
        catch
        end
        fprintf('----------------------------------------------------------\n')
        fprintf('Data saved to file. Have a nice day!\n')
    end
    
    function[E] = gen_e
        % Generates the value E for each connection type as specified on page 9 of
        % Humphries, Wood & Gurney (2010)

        % Set maximum distance to greater than striatal edge length to allow for 
        % corner-to-corner checks
        d = 1:2000;

        % Values from page 10, table 4
        alpha = [0.511, -0.921, -0.695, 1.322];
        beta  = [1.033, 1.033, 1.38, 2.4];
        gamma = [0.042, 0.042, 0.057, 0.016];
        delta = [26.8, 26.8, 15.6, 43.3];
        eta   = [0.0039, 0.0039, 0.0036, 0.0029];

        % Page 9, equation 20
        for c = 1:numel(alpha) 
            E(c,:) = exp(-alpha(c) - beta(c) .* (1 - exp(-gamma(c) .* (d - delta(c)))) .* exp(eta(c) .* d));
        end
    end
end

