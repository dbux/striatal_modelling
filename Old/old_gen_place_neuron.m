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