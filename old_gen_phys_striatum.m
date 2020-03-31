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
    fprintf('\nInitializing ?m^3 striatum with %d MSNs + %d%% FSIs? ', attr.size, pop.msn, attr.fsi_pct)
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
    fprintf('Saving striatum data? ')
end
striatum.dirname = [attr.path num2str(strinf.id) '_' ...
    num2str(pop.msn) '+' num2str(pop.fsi) '_' num2str(attr.ch_all) 'CH_' num2str(flags.phys_ch) 'sep_' num2str(attr.ch_overlap) 'overlap/'];
mkdir(striatum.dirname);
filename = [striatum.dirname '/striatum.mat'];
save(filename, 'striatum', 'strinf', 'attr');
if flags.progress
    fprintf('took %1.2f seconds. Done!\n', toc(timer.save))
end