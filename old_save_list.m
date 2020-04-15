function[] = save_list(listpath, listfile, name, flags)
% Given a pathname, a connection list and name, this will save the
% connection list to a CSV and optionally convert it to binary format for
% direct import to SpineCreator

% Save connection list to CSV
filename = ['conn_', name.src, '_to_', name.dst, '_', name.syn];
fid = fopen([listpath, filename, '.csv'], 'w');

% With or without connection delays as appropriate
if size(listfile,2) == 3
    fprintf(fid, '%d, %d, %1.1f\r\n', transpose(listfile));
elseif size(listfile,2) == 2
    fprintf(fid, '%d, %d\r\n', transpose(listfile));
else
    error('Incorrect number of columns in connection %s', filename);
end
fclose(fid);

% Export connection list as binary file if need
if flags.binary
    bin_convert(listfile, strcat(listpath, filename));
end
