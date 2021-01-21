function[striatum] = gen_phys_striatum(attr)
    % Generates a model striatum with topology based on the description in
    % page 7 of Humphries, Wood & Gurney (2010)
    
    % Extract necessary attributes
    phys = attr.phys;
    conn = attr.conn;
    flags = attr.flags;

    pop.msn = floor(round(phys.size^3 * phys.msn_density) / 1e9);   % Number of MSNs to place
    while mod((pop.msn/2), conn.ch_all) ~= 0                        % Ensure the number of MSNs will fit nicely into the number of channels
        pop.msn = pop.msn + 1;
    end

    pop.max = round(pop.msn + (pop.msn/10));                        % Maximum possible population (for matrix preallocation only)
    strinf.centre = phys.size/2;                                    % Centre point of striatum

    % How neurons will be represented numerically
    msn = 1;                
    fsi = 3;

    % Preallocate striatal array - initialize as uint8 to save memory and time
    % Striatal centre arrays small enough to not require preallocation
    striatum.main = zeros(phys.size, phys.size, phys.size, 'uint16');
    striatum.linear = zeros(pop.max,1, 'uint8');
    striatum.linear_centre = [];
    striatum.neurons = zeros(pop.max,3, 'uint16');
    striatum.neurons_centre = [];
    strinf.ptr = 1;                                                 % Points to the next available striatum entry

    if flags.progress
        timer.str = tic;
        reverseStr = '';                                            % Required string for progress indicator
        fprintf('Initializing %dμm³ striatum with %d MSNs + %d%% FSIs… ', phys.size, pop.msn, phys.fsi_pct)
    end

    % Place each neuron individually
    for i = 1:pop.msn  
        [striatum, strinf] = place_neuron(striatum, strinf, phys, msn);

        % Roll the dice to see if we will also be generating an FSI
        % FSI creation procedure is the same as above
        if randi(100) <= phys.fsi_pct
            [striatum, strinf] = place_neuron(striatum, strinf, phys, fsi);       
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
end

%% FUNCTIONS
function[striatum, strinf] = place_neuron(striatum, strinf, phys, n_type)
    % Given a striatum structure and associated information, will place a
    % neuron in the physical striatum ensuring a minimum distance from all
    % other neurons and return the co-ordinates of that neuron

    % !IMPORTANT! This program requires the 'Geom3D' toolbox available from:
    % http://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
    % Specifically, the files 'distancePoints3d.m' and 'vectorNorm3d.m' must be
    % available for calculating neuron distances.

    % Reset 'neuron OK' flag before entering while-loop
    n_ok = 0;

    while n_ok == 0
        % Generate random xyz co-ordinates for neuron placement
        n_loc = randi(phys.size, 1, 3);

        % Assume that the co-ords are OK unless we decide otherwise
        n_ok = 1;

        % If those co-ords are already occupied, generate a new set
        while striatum.main(n_loc(1), n_loc(2), n_loc(3)) ~= 0
            n_loc = randi(phys.size, 1, 3);
        end

        % Create a set of 'buffer' co-ordinates that extend the minimum
        % separation distance away from the potential neuron co-ords in all directions
        for j = 1:3
            if n_loc(j) - phys.min_dist < 1
                buffer{j}(1) = 1;
            else
                buffer{j}(1) = n_loc(j) - phys.min_dist;
            end

            if n_loc(j) + phys.min_dist > phys.size
                buffer{j}(2) = phys.size;
            else
                buffer{j}(2) = n_loc(j) + phys.min_dist;
            end
        end

        % Test to see if any neuron already in this area would be too close
        % to the new neuron. Test is for distance, not just presence,
        % because neuron may be close in only one or two dimensions

        % Define buffer zone
        n_buffer = striatum.main(...
                buffer{1}(1):buffer{1}(2), ...
                buffer{2}(1):buffer{2}(2), ...
                buffer{3}(1):buffer{3}(2));

        % If any neurons are present in the buffer zone
        if any(n_buffer, 'all')
            % Get co-ordinates of all neurons in buffer zone
            [x, y, z] = ind2sub(size(n_buffer), find(n_buffer));
            for n = 1:numel(x)               
                % If minimum distance is violated then set 'neuron not OK'
                if distancePoints3d(n_loc, [x(n), y(n), z(n)]) < phys.min_dist         
                    n_ok = 0;
                end
            end
        end

    end

    % If the new neuron is close to the centre of the striatum, add it to the list of central neurons
    if distancePoints3d([strinf.centre, strinf.centre, strinf.centre], n_loc) <= phys.centre_rad
        striatum.linear_centre  = [striatum.linear_centre ; strinf.ptr, n_type];
        striatum.neurons_centre = [striatum.neurons_centre ; n_loc];
    end

    % Everything checks out! Add the neuron to the striatum and the neuron location register
    striatum.main(n_loc) = n_type;
    striatum.linear(strinf.ptr) = n_type;
    striatum.neurons(strinf.ptr,:) = n_loc;
    strinf.ptr = strinf.ptr + 1;
end