function[] = bin_convert(data, filename)
% Given a two- or three-column array of data representing input neuron,
% output neuron and (optionally) delay, this function converts the data to 
% binary format so that it can be dropped straight into the SpineCreator
% model directory.

intArray = int32([data(1:end,1)' ; data(1:end,2)']');
if size(data, 2) == 3
    floatArray = single(data(1:end,3));
end

fileID = fopen(strcat(filename,'.bin'),'w');

if exist('floatArray', 'var')
    for i = 1:size(data,1)
        fwrite(fileID, intArray(i,:), 'int32');
        fwrite(fileID, floatArray(i), 'single');
    end
else
    for i = 1:size(data,1)
        fwrite(fileID, intArray(i,:), 'int32');
    end
end
fclose(fileID);