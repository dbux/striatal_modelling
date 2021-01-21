function[] = save_list(path, file, name, flags)
    % Given a pathname, a connection list and name, this will save the
    % connection list to a CSV and optionally convert it to binary format for
    % direct import to SpineCreator

    % Save connection list to CSV
    filename = ['conn_', name.src, '_to_', name.dst, '_', name.syn];
    fid = fopen(fullfile(path, [filename, '.csv']), 'w');

    % With or without connection delays as appropriate
    if size(file,2) == 3
        fprintf(fid, '%d, %d, %1.1f\r\n', transpose(file));
    elseif size(file,2) == 2
        fprintf(fid, '%d, %d\r\n', transpose(file));
    else
        error('Incorrect number of columns in connection %s', filename);
    end
    fclose(fid);

    % Export connection list as binary file if needed
    if flags.binary
        bin_convert(file, fullfile(path, filename));
    end
    
    
    function[] = bin_convert(data, filename)
        % Given a two- or three-column array of data representing input neuron,
        % output neuron and (optionally) delay, this function converts the data to 
        % binary format so that it can be dropped straight into the SpineCreator model directory.

        intArray = int32([data(1 : end, 1)' ; data(1 : end, 2)']');
        if size(data, 2) == 3
            floatArray = single(data(1 : end, 3));
        end

        fileID = fopen([filename, '.bin'], 'w');

        if exist('floatArray', 'var')
            for i = 1:size(data,1)
                fwrite(fileID, intArray(i, :), 'int32');
                fwrite(fileID, floatArray(i), 'single');
            end
        else
            for i = 1:size(data,1)
                fwrite(fileID, intArray(i, :), 'int32');
            end
        end
        
        fclose(fileID);
    end       
end


