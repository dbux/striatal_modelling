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
%
% Physical striatum attributes
attr.size        = 600;    % Size of model striatum each side (μm) (Was 250μm in Humphries et al. 2009)
attr.min_dist    = 10;      % Minimum distance between neurons (μm)
attr.centre_rad  = 75;      % Radius of central region free from edge effects (μm)
attr.msn_density = 84900;   % Number of MSNs to place per mm^3 (should be 84,900)
attr.fsi_pct     = 1;       % Percentage of MSNs to be added as FSIs  

% Connectivity attributes
attr.ch_all     = 2;        % Total number of channels (not including background)
attr.ch_overlap = 0;        % Overlap (in % of striatal size) for two channels in physical model
attr.bkg_msn    = 0;        % Percentage of MSNs to receive only background noise. Leave at 0 for no background.
attr.bkg_fsi    = 0;        % Percentage of FSIs to receive only background noise. Leave at 0 for no background.

% Process flags
flags.phys_ch   = 1;        % Separate channels based on physical location? (TWO CHANNELS ONLY)
flags.debug     = 0;        % Show detailed information during initialization?
flags.progress  = 0;        % Show progress indicator? (Set to 0 for Iceberg)
flags.save      = 1;        % Save connection lists to disk?
flags.binary    = 1;        % Save binary versions of connection lists?
flags.density   = 1;        % Create input lists for varying neural densities?
flags.overlap   = 1;        % Create inputs lists for varying channel overlap?

% % % % % % % % % %
% BEGIN
%
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

    % Create corticostriatal connections and output all connection lists
    % TODO: Loop this cleanly
    for i = 0:10:90
        for j = 0:10:90
            % Iterate active neuron density using flags.density and flags.overlap
            attr.bkg_msn = i;
            attr.bkg_fsi = j;

            [connections, list] = gen_phys_connlist(striatum, connections, attr, flags);
        end
    end
    
    for k = -90:10:90
        % Iterate channel overlap
        attr.ch_overlap = k;

        [connections, list] = gen_phys_connlist(striatum, connections, attr, flags);
    end
end

% % % % % % % % % %
% FUNCTIONS
%
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
        num2str(pop.msn) '+' num2str(pop.fsi) '_' num2str(attr.ch_all) 'CH/'];
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
    %     fprintf('\nEstimated completion time: %1.2f - %1.2f minutes… ', (0.09*length(striatum.linear))/60, (0.15*length(striatum.linear))/60)
    end

    % For each neuron...
    for i = 1:length(striatum.linear)
        timer.neuron = tic;

        % Connection information
        coninf.count    = 0;    % Counts number of connections created 
        coninf.dist_all = [];	% Distance from this neuron to all other neurons
        coninf.exp_all  = [];   % Expected number of connections of all types to all other neurons 
        coninf.exp_this = [];   % Expected number of connections of appropriate types to all other neurons
        coninf.trials   = [];   % Number of connection trials to perform for all potential connections
        coninf.prob     = [];   % Probability of making a connection at each trial to all other neurons

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
    fprintf(fc, '\nConnection type, No. contacts, STD, Distance (µm), STD');

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

function[connections, list] = gen_phys_connlist(striatum, connections, attr, flags)
    % Assign MSNs as D1 or D2 and to a particular channel

    % Append directory for connection lists
    listpath = [striatum.dirname 'connection_lists/'];
    mkdir(listpath);

    % Start list timer
    timer.list = tic;

    % Create a list of all MSN, FSI, and D1 / D2 neuron IDs (without connections)
    % Column 1 is the unique neuron ID in the striatum generated by MatLab
    % Column 2 is the neuron's ID within that group required by SpineCreator
    list.msn = unique(connections.msnmsn(:,1));
    list.d1 = [list.msn(1 : ceil(length(list.msn) / 2))' ; 0 : ceil(length(list.msn) / 2) - 1]';
    list.d2 = [list.msn(length(list.d1) + 1 : end)' ; 0 : floor(length(list.msn) / 2) - 1]';
    list.fsi = sort(unique([unique(connections.fsifsi(:,1)) ; unique(connections.fsimsn(:,1))]));
    list.fsi(:,2) = 0:length(list.fsi) - 1;

    % Separate all MSN_MSN connections into D_type lists
    timer.chans = tic;
    if flags.progress
        fprintf('\nSeparating connections into channels… ');
    end
    connections.d1.all = connections.msnmsn((connections.msnmsn(:, 1) <= max(list.d1(:, 1))), :);
    connections.d1.d1  = connections.d1.all((connections.d1.all(:, 2) <= max(list.d1(:, 1))), :);
    connections.d1.d2  = connections.d1.all((connections.d1.all(:, 2) > max(list.d1(:, 1))), :);

    connections.d2.all = connections.msnmsn((connections.msnmsn(:, 1) > max(list.d1(:, 1))), :);
    connections.d2.d1  = connections.d2.all((connections.d2.all(:, 2) <= max(list.d1(:, 1))), :);
    connections.d2.d2  = connections.d2.all((connections.d2.all(:, 2) > max(list.d1(:, 1))), :);

    connections.fsi.d1 = connections.fsimsn((connections.fsimsn(:, 2) <= max(list.d1(:, 1))), :);
    connections.fsi.d2 = connections.fsimsn((connections.fsimsn(:, 2) > max(list.d1(:, 1))), :);
    if flags.progress
        fprintf('done! (%3.2fs)\n', toc(timer.chans))
    end

    % It's useful to know how many of each neuron type there are
    num.d1 = size(list.d1, 1);
    num.d2 = size(list.d2, 1);
    num.msn = size(list.msn, 1);
    num.fsi = size(list.fsi, 1);
    num.gap = size(connections.gap, 1);

    % (Approximate) number of MSNs of each type to leave as background only
    num.bkg_msn = round(num.msn * (attr.bkg_msn / 100) / 2);

    % (Approximate) number of FSIs to leave as background only
    num.bkg_fsi = round(num.fsi * (attr.bkg_fsi / 100));

    % Number of MSNs of each type and FSIs to put in each channel
    % num.msn_ch = floor((num.msn / 2 - num.bkg) / attr.ch_all);
    num.d1_ch = (num.d1 - num.bkg_msn) / attr.ch_all;
    num.d2_ch = (num.d2 - num.bkg_msn) / attr.ch_all;
    num.fsi_ch = num.fsi - num.bkg_fsi;

    %% CORTICO-STRIATAL channel connections - Uses SpineCreator neuron IDs
    if flags.progress
        fprintf('\nCreating connection lists:\n');
        fprintf('1) Cortical channel connections… ')
    end
    timer.conn1 = tic;

    % Cortico-striatal connections differ based on number of input channels
    if attr.ch_all == 1
        % In the single-channel model it doesn't matter which neurons don't
        % receive connections since neuron ID is not associated with location

        connections.cortex.ch1.d1  = [0 : num.d1_ch - 1; 0 : num.d1_ch - 1]';
        connections.cortex.ch1.d2  = [(0 : num.d2_ch - 1) + num.d1_ch ; 0 : num.d2_ch - 1]';
        connections.cortex.ch1.fsi = [0 : num.fsi_ch - 1; 0 : num.fsi_ch - 1]';

        % Save connection lists
        if flags.save
            % For both D1 and D2 MSNs
            for j = 1:2           
                d_dst = sprintf('d%d', j);
                name.src = sprintf('CH1_input');
                name.dst = sprintf('Striatum_D%d', j);

                % Create both AMPA and NMDA connections to MSNs
                for k = 0:1          
                    % Channel input

                    name.syn = sprintf('syn%d_bkMSN%d_bkFSI%d', k, attr.bkg_msn, attr.bkg_fsi);
                    save_list(listpath, connections.cortex.ch1.(d_dst), name, flags);
                end
            end

            % FSI populations only use a single synapse
            name.dst = 'Striatum_FSI';
            name.syn = sprintf('syn0_bkMSN%d_bkFSI%d', attr.bkg_msn, attr.bkg_fsi);
            save_list(listpath, connections.cortex.ch1.fsi, name, flags);
        end

    elseif attr.ch_all == 2

        % TODO: EXPORT CH LISTS AS CO_ORDS FOR PLOTTING

        % Physically partition a two-channel striatum
        if flags.phys_ch  
            for i = 1:2
                % Set dynamic structure fieldname
                msn = sprintf('d%d', i);

                % MSNs assigned to channel 1 or 2 based on physical location on
                % striatal X-axis, modified based on overlap value. Negative
                % overlap creates a background-only gap between channels,
                % positive overlap creates a region with MSNs in both channels
%                 list.ch1.(msn) = list.(msn)(striatum.neurons(list.(msn)(:,1), 1) <= attr.size / 2 + (attr.size * attr.ch_overlap) / 2, :);
%                 list.ch2.(msn) = list.(msn)(striatum.neurons(list.(msn)(:,1), 1) >  attr.size / 2 + (attr.size * attr.ch_overlap) / 2, :);     
                list.ch1.(msn) = list.(msn)(striatum.neurons(list.(msn)(:,1), 1) <= (attr.size / 2) + (attr.size * (attr.ch_overlap / 100) / 2), :);
                list.ch2.(msn) = list.(msn)(striatum.neurons(list.(msn)(:,1), 1) >= (attr.size / 2) - (attr.size * (attr.ch_overlap / 100) / 2), :); 

                for j = 1:attr.ch_all
                    % Set dynamic structure fieldname
                    ch = sprintf('ch%d', j);

                    % Trim MSN and FSI lists according to requested background-only percentage
                    list.(ch).(msn) = list.(ch).(msn)(1:end - (floor(size(list.(ch).(msn), 1) * (attr.bkg_msn / 100))),:);
                    list.(ch).fsi = list.fsi(1:end - (floor(size(list.fsi, 1) * (attr.bkg_fsi / 100))),:);

                    % Create channel connections to striatum
                    % FROM: Each cortical channel (-1 for SpineCreator 0-indexing)
                    % TO:   Each MSN or FSI in each channel
                    % DELAY:N/A
%                     connections.cortex.(ch).(msn) = gen_list_allall(attr.ch_inputs, list.(ch).(msn)(:,2));
                    connections.cortex.(ch).(msn) = [0 : length(list.(ch).(msn)) - 1 ; list.(ch).(msn)(:,2)']';
%                     connections.cortex.(ch).fsi = gen_list_allall(attr.ch_inputs, list.(ch).fsi(:,2));
                    connections.cortex.(ch).fsi = [0 : length(list.(ch).fsi) - 1 ; list.(ch).fsi(:,2)']';
                end
            end

    %     % Normal procedure for non-partitioned multi-channel striatum    
        else  
            fprintf('\nNew connections for non-partioned striatum not yet done!')
    %         % Position of the first and last neuron in the channel
    %         ch_first = 1;
    %         ch_last = num.msn_ch;
    % 
    %         for i = 1:attr.ch_all
    %             % Set dynamic structure fieldname
    %             ch = sprintf('ch%d', i);
    % 
    %             % The neurons structure contains a list of members of each channel
    %             list.(ch).d1 = list.d1(ch_first : ch_last, :);
    %             list.(ch).d2 = list.d2(ch_first : ch_last, :);
    % 
    %             % Create channel connections to striatum
    %             % FROM: Each cortical channel (-1 for SpineCreator 0-indexing)
    %             % TO:   Each MSN or FSI in each channel
    %             % DELAY:N/A
    %             connections.cortex.(ch).d1 = [0 : num.msn_ch - 1 ; list.(ch).d1(:,2)']';
    %             connections.cortex.(ch).d2 = [0 : num.msn_ch - 1 ; list.(ch).d2(:,2)']';
    %             connections.cortex.(ch).fsi= [0 : num.fsi - 1    ; list.fsi(:,2)'    ]'; 
    % 
    %             % Create channel connections to basal ganglia
    %             % FROM: Each cortical channel (-1 for SpineCreator 0-indexing)
    %             % TO:   Single BG neuron in each channel
    %             % DELAY:N/A
    %             if num.msn_ch > attr.max_bg
    %                 connections.cortex.(ch).bg = [0 : attr.max_bg - 1 ; zeros(1, attr.max_bg) + (i - 1)]'; 
    %             else
    %                 connections.cortex.(ch).bg = [0 : num.msn_ch - 1 ; zeros(1, num.msn_ch) + (i - 1)]';  
    %             end
    % 
    %             % Increment channel start/end markers for the next channel
    %             ch_first = ch_last + 1;
    %             ch_last = ch_last + num.msn_ch;
    %         end
        end

        % Save connection lists
        if flags.save
            for i = 1:attr.ch_all
                % Set dynamic structure fieldname
                ch = sprintf('ch%d', i);       
                name.src = sprintf('CH%d_input', i);

                % For both D1 and D2 MSNs
                for j = 1:2           
                    d_dst = sprintf('d%d', j);
                    name.dst = sprintf('Striatum_D%d', j);

                    % Create both AMPA and NMDA connections to MSNs
                    for k = 0:1
                        name.syn = sprintf('syn%d_bkMSN%d_bkFSI%d_ovlp%d', k, attr.bkg_msn, attr.bkg_fsi, attr.ch_overlap);                
                        save_list(listpath, connections.cortex.(ch).(d_dst), name, flags);
                    end
                end

                % FSIs only use a single synapse
                name.dst = 'Striatum_FSI';
                name.syn = sprintf('syn0_bkMSN%d_bkFSI%d_ovlp%d', attr.bkg_msn, attr.bkg_fsi, attr.ch_overlap);
                save_list(listpath, connections.cortex.(ch).fsi, name, flags);
            end
        end

    elseif attr.ch_all > 2
        fprintf('3+ channel connectivity not yet done')
    end


    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn1))
    end


    %% STRIATAL GABA connections - Converts MatLab to SpineCreator neuron IDs
    if flags.progress
        fprintf('2) Striatal GABA connections… ')
    end
    timer.conn3 = tic;

    % From both D1 and D2 MSNs
    for i = 1:2
        d_src = sprintf('d%d', i); 
        name.src = sprintf('Striatum_D%d', i);

        % GABA projections always use syn0
        name.syn = 'syn0';

        % To both D1 and D2 MSNs
        for j = 1:2 
            d_dst = sprintf('d%d',j);
            name.dst = sprintf('Striatum_D%d', j);

            % If connections already exist, don't recreate
            if ~exist(sprintf([listpath, 'conn_', name.src, '_to_', name.dst, '_', name.syn, '.csv']), 'file')        
                % Convert MatLab neuron IDs to SpineCreator IDs
                [~, src] = ismember(connections.(d_src).(d_dst)(:,1), list.(d_src)(:,1));
                [~, dst] = ismember(connections.(d_src).(d_dst)(:,2), list.(d_dst)(:,1));

                % Create list of undirected MSN-MSN connections
                % FROM: D1 or D2 MSNs
                % TO:   D1 or D2 MSNs
                % DELAY:As defined in connections.(d_src).(d_dst)(:,3)
                connections.gaba.(d_src).(d_dst) = ...
                    [list.(d_src)(src,2)' ; list.(d_dst)(dst,2)' ; connections.(d_src).(d_dst)(:,3)']';

                % Save connection lists
                if flags.save
                    save_list(listpath, connections.gaba.(d_src).(d_dst), name, flags);
                end
            end
        end

        % From FSIs to MSNs
        f_dst = sprintf('d%d', i); 
        name.src = 'Striatum_FSI';
        name.dst = sprintf('Striatum_D%d', i);

        % If connections already exist, don't recreate
        if ~exist(sprintf([listpath, 'conn_', name.src, '_to_', name.dst, '_', name.syn, '.csv']), 'file') 
            % Convert MatLab neuron IDs to SpineCreator IDs
            [~, src] = ismember(connections.fsi.(f_dst)(:,1), list.fsi(:,1));
            [~, dst] = ismember(connections.fsi.(f_dst)(:,2), list.(f_dst)(:,1));

            % Create list of undirected FSI-MSN connections
            % FROM: FSIs
            % TO:   D1 or D2 MSNs
            % DELAY:As defined in connections.fsi.(f_dst)(:,3)
            connections.gaba.fsi.(f_dst) = [list.fsi(src,2)' ; list.(f_dst)(dst,2)' ; connections.fsi.(f_dst)(:,3)']';

            % Save connection lists
            if flags.save
                save_list(listpath, connections.gaba.fsi.(f_dst), name, flags);
            end   
        end
    end

    % From FSIs to FSIs (GABA)
    name.src = 'Striatum_FSI';
    name.dst = 'Striatum_FSI';
    name.syn = 'syn0'; 

    % If connections already exist, don't recreate
    if ~exist(sprintf([listpath, 'conn_', name.src, '_to_', name.dst, '_', name.syn, '.csv']), 'file') 
        % Convert MatLab neuron IDs to SpineCreator IDs
        [~, src] = ismember(connections.fsifsi(:,1), list.fsi(:,1));
        [~, dst] = ismember(connections.fsifsi(:,2), list.fsi(:,1));

        % Create list of undirected FSI-FSI GABA connections
        % FROM: FSIs
        % TO:   FSIs
        % DELAY:As defined in connections.fsifsi(:,3)
        connections.gaba.fsi.fsi = [list.fsi(src,2)' ; list.fsi(dst,2)' ; connections.fsifsi(:,3)']';

        % Save connection lists
        if flags.save
            save_list(listpath, connections.gaba.fsi.fsi, name, flags);
        end  

        % From FSIs to FSIs (Gap)
        % Convert MatLab neuron IDs to SpineCreator IDs
        try
            [~, src] = ismember(connections.gap(:,1), list.fsi(:,1));
            [~, dst] = ismember(connections.gap(:,2), list.fsi(:,1));
        catch
        end

        % Create list of FSI-FSI gap connections
        % Gap junctions have a nonstandard format
        % FROM: FSIs
        % TO:   FSIs
        % DELAY:N/A
        try
            connections.gap_sc.in1 =  [list.fsi(src,2)' ; 0 : num.gap - 1]';
            connections.gap_sc.in2 =  [list.fsi(dst,2)' ; 0 : num.gap - 1]';
            connections.gap_sc.out1 = [0 : num.gap - 1  ; list.fsi(src,2)']';
            connections.gap_sc.out2 = [0 : num.gap - 1  ; list.fsi(dst,2)']';        
        catch
            error('Could not create gap junctions!')
        end

        % Save connection lists
        if flags.save
            for i = 0:1
                g_in = sprintf('in%d', i + 1); 
                g_out = sprintf('out%d', i + 1); 

                name.syn = sprintf('syn%d', i);

                name.src = 'Striatum_FSI';               
                name.dst = 'FSI_GAP';
                save_list(listpath, connections.gap_sc.(g_in), name, flags);

                name.src = 'FSI_GAP';
                name.dst = 'Striatum_FSI';
                save_list(listpath, connections.gap_sc.(g_out), name, flags);        
            end
        end
    end

    if flags.progress
        fprintf('done! (%1.2fs)\n', toc(timer.conn3))
    end

    %% Save neuron list and connections to disk
    timer.save = tic;
    if flags.progress
        fprintf('\nSaving neuron and connection data… ')
    end

    filename = [striatum.dirname '/connections.mat'];
    save(filename, 'connections', '-v7.3');

    filename = [striatum.dirname '/list.mat'];
    save(filename, 'list');
    
    % Get physical co-ordinates of all MSNs and FSIs
    striatum.msn.d1  = striatum.neurons(list.d1(:, 1), :);
    striatum.msn.d2  = striatum.neurons(list.d2(:, 1), :);
    striatum.msn.all = [striatum.msn.d1 ; striatum.msn.d2];
    striatum.fsi.all = striatum.neurons(list.fsi(:, 1), :);
    
    striatum.fsi.active    = [];
    striatum.fsi.inactive  = [];
    striatum.msn.active.d1 = [];
    striatum.msn.active.d2 = [];
    
    for i = 1:attr.ch_all
        % Set dynamic structure fieldname
        ch = sprintf('ch%d', i);
        
        % Get active FSIs
        [~, idf] = ismember(connections.cortex.(ch).fsi(:, 2), list.fsi(:, 2));
        striatum.fsi.active = unique([striatum.fsi.active ; striatum.fsi.all(idf, :)], 'rows');
        
        striatum.msn.(ch).all = [];

        % Get active MSNs
        for j = 1:2   
            % Set dynamic structure fieldname
            dx = sprintf('d%d', j);
            
            [~, idm] = ismember(connections.cortex.(ch).(dx)(:, 2), list.(dx)(:, 2));
            striatum.msn.(ch).(dx) = striatum.neurons(list.(dx)(idm, 1), :);
            striatum.msn.(ch).all  = [striatum.msn.(ch).all ; striatum.msn.(ch).(dx)];
            striatum.msn.active.(dx) = [striatum.msn.active.(dx) ; striatum.msn.(ch).(dx)];
        end
    end

    % Get inactive neurons
    striatum.msn.inactive.d1  = striatum.msn.d1(~ismember(striatum.msn.d1, striatum.msn.active.d1, 'rows'), :);
    striatum.msn.inactive.d2  = striatum.msn.d2(~ismember(striatum.msn.d2, striatum.msn.active.d2, 'rows'), :);
    striatum.msn.inactive.all = [striatum.msn.inactive.d1 ; striatum.msn.inactive.d2];
    striatum.fsi.inactive     = striatum.fsi.all(~ismember(striatum.fsi.all, striatum.fsi.active, 'rows'), :);

    % Export co-ordinates of all MSNs and FSIs
    export_striatum(striatum.msn, 'striatum.msn', striatum.dirname)
    export_striatum(striatum.fsi, 'striatum.fsi', striatum.dirname)
    
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
                cell2csv([s_dir, '/', curr_name, '.csv'], ...
                    [{[curr_name, '_X'], [curr_name, '_Y'], [curr_name, '_Z']} ; ...
                    num2cell(s_field)]);
            end
        end   
    end
        
    if flags.progress
        fprintf('took %1.2f minutes. All done!\n', toc(timer.save)/60)
    end

    if flags.progress
        fprintf('All connection lists generated in %1.2f minutes.\n', toc(timer.list) / 60)
    end
end

