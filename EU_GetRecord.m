% Read data from the EU database files. Use metadata from http://epilepsiae.uniklinik-freiburg.de
% to annotate these files with the clinical and subclincal seizure onset and end times.
% GOL 2017 - University of Toronto
function   [ record ] = EU_GetRecord(EU_data_path, patient_id, record_idx, channels, downsample_rate)

    %% Load paths for patient files use .mat file to save time after first run
    [current_dir,name,ext] = fileparts(mfilename('fullpath'));
    metadata_dir = [current_dir, '/metadata/'];  
    
    fname = [metadata_dir, patient_id, '_record.mat'];
    if exist(fname, 'file') == 0 
        [record_list, seizure_records, onset_samples, offset_samples, sc_records, sc_onset_samples, sc_offset_samples, ieeg_onset, ieeg_offset]= EU_GetMetadata(EU_data_path, patient_id);

        % Write a mat file with records info
        save(fname, ...
            'record_list',...
            'ieeg_onset',...
            'ieeg_offset', ...                  
            'seizure_records',...
            'onset_samples',...
            'offset_samples', ...
            'sc_records', ...
            'sc_onset_samples', ...
            'sc_offset_samples', ...
            '-v7.3' ); % Use when filename is stored in a variable
    else
        load(fname);
    end

    %% Import iEEG data
    selected_record = record_list{record_idx};
    file=bin_file(selected_record);

    % need to use def_data_access before get_bin_signals, otherwise
    % a_channs_cell is not set, and this fails in get_bin_signals
    %             if isempty(self.a_file_elec_cell) || isempty(self.a_channs_cell)
    %                 return;
    %             end
    %def_data_access(self, wsize, step, channs_cell, offset)
    % however, this function can just set it as:  self.a_channs_cell = self.a_file_elec_cell;
    file.a_channs_cell = file.a_file_elec_cell;
    
    chan_names = file.a_file_elec_cell;
    %chan_names = cellstr(chan_names);
    
    %% Convert channel names into the corresponding index column
    n_select_chans = size(channels,1);
    n_avail_chans = size(chan_names,2);
    
    if(strcmp(channels(1),'ALL'))
        chan_idxs = 1:n_avail_chans;
        n_select_chans = n_avail_chans;
        channels = chan_names;
    else
        chan_idxs = zeros(n_select_chans,1);
        for c = 1:n_select_chans
            index = find(strcmp(chan_names, channels(c,:)));
            chan_idxs(c) = index;
        end
    end
   
    %% Get data
    fprintf('Loading data...');
    
    tmp_dir = ['./tmp_' , patient_id];
    mkdir(tmp_dir);
    workingDir = pwd; 
    cd(tmp_dir)    % to avoid overlapping temp files
    n_samps = file.a_n_samples;
    t = tic;
    ts_data=file.get_bin_signals(1, n_samps );
    elapsed = toc(t);
    fprintf(' (time: %f seconds)\n', elapsed);
    cd(workingDir)
  
    signal = ts_data(chan_idxs,:)';
    
    sig_ds = zeros( ceil(size(signal,1)/downsample_rate), ceil(size(signal,2)) );
    if (downsample_rate > 1)
            for c = 1:n_select_chans
                sig_ds(:,c) = decimate(signal(:,c),downsample_rate);
            end
            signal = sig_ds;
    end
    
    Fs = file.a_samp_freq / downsample_rate;
    
    if (Fs > 256)
        fprintf('Warning: Sample rate is not 256 - is this correct? (enter to continue)');
        pause;
    end

    
    %% Check if a seizure record
    sz_idx = find(strcmp(seizure_records,selected_record));
    %sz_idx = find(contains());
    
    if (sz_idx)
        for i = 1:length(sz_idx)
        sample_sz_onset{i} = onset_samples{sz_idx(i)} / downsample_rate;
        sample_sz_end{i} = offset_samples{sz_idx(i)} / downsample_rate;           
        end        
    else
        sample_sz_onset = {0};
        sample_sz_end = {0};
    end
    
    %% Check if a subclinical seizure record
    scsz_idx = find(strcmp(sc_records,selected_record));
    %sz_idx = find(contains());
    
    if (scsz_idx)
        for i = 1:length(scsz_idx)
            sample_scsz_onset{i} = sc_onset_samples{scsz_idx(i)} / downsample_rate;
            sample_scsz_end{i} = sc_offset_samples{scsz_idx(i)} / downsample_rate;       
        end
    else
        sample_scsz_onset = {0};
        sample_scsz_end = {0};
    end  
    

    record.signal = signal;
    record.record_idx = record_idx;
    record.Fs = Fs;
    record.patient_id = patient_id;
    record.channels = channels;
    record.selected_record = selected_record;
    record.sample_sz_onset = sample_sz_onset;
    record.sample_sz_end = sample_sz_end;
    record.sample_scsz_onset =  sample_scsz_onset;
    record.sample_scsz_end = sample_scsz_end;
    record.isSZ = sz_idx;
    record.isSCSZ = scsz_idx;
    record.ieeg_onset = ieeg_onset{record_idx};
    record.ieeg_offset = ieeg_offset{record_idx};
end