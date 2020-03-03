% load('/media/dbuxton/Data/Dropbox/University/list.mat')
% load('/media/dbuxton/Data/Dropbox/University/striatum.mat')

listpath = './connection_lists/';
mkdir(listpath);

num.msn = sum(striatum.linear(:) == 1);
num.fsi = sum(striatum.linear(:) == 3);

% It's useful to know how many of each neuron type there are
num.d1 = ceil(num.msn / 2);
num.d2 = floor(num.msn / 2);

for i = 0:10:90
    for j = 0:10:90
        % Iterate active neuron density
        attr.bkg_msn = i;
        attr.bkg_fsi = j;
        
        % (Approximate) number of MSNs of each type to leave as background only
        num.bkg_msn = round(num.msn * (attr.bkg_msn / 100) / 2);

        % (Approximate) number of FSIs to leave as background only
        num.bkg_fsi = round(num.fsi * (attr.bkg_fsi / 100));

        % Number of MSNs of each type and FSIs to put in each channel
        % num.msn_ch = floor((num.msn / 2 - num.bkg) / attr.ch_all);
        num.d1_ch = (num.d1 - num.bkg_msn) / attr.ch_all;
        num.d2_ch = (num.d2 - num.bkg_msn) / attr.ch_all;
        num.fsi_ch = num.fsi - num.bkg_fsi;

        connections.cortex.ch1.d1 = [0 : num.d1_ch - 1 ; 0 : num.d1_ch - 1]';
        connections.cortex.ch1.d2 = [(0 : num.d2_ch - 1) + ceil(num.msn / 2); 0 : num.d2_ch - 1]';
        connections.cortex.ch1.fsi = [0 : num.fsi_ch - 1 ; 0 : num.fsi_ch - 1]';

        % For both D1 and D2 MSNs
        for k = 1:2           
            d_dst = sprintf('d%d', k);
            name.src = sprintf('CH1_input');
            name.dst = sprintf('Striatum_D%d', k);

            % Create both AMPA and NMDA connections to MSNs
            for l = 0:1          
                % Channel input

                name.syn = sprintf('syn%d_bkMSN%d_bkFSI%d', l, attr.bkg_msn, attr.bkg_fsi);
                save_list(listpath, connections.cortex.ch1.(d_dst), name, flags);
            end
        end

        % FSI populations only use a single synapse
        name.dst = 'Striatum_FSI';
        name.syn = sprintf('syn0_bkMSN%d_bkFSI%d', attr.bkg_msn, attr.bkg_fsi);
        save_list(listpath, connections.cortex.ch1.fsi, name, flags);
    end
end