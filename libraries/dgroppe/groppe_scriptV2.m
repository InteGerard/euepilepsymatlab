%% Loop over header files and order them from start to finish
%patient_root='/Volumes/Seagate Expansion Drive/EU/inv/pat_FR_253/adm_253102/';

patient_name = 'pat_FR_253';

patient_root=['/media/gerard/My Book/EU/inv/', patient_name]; %pat_FR_253/adm_253102'; % Edit this for each new patient or drive location
%patient_root=['/media/gerard/My Book/EU/inv/pat_FR_253/adm_253102']; % Edit this for each new patient or drive location

% Find adm_ directory and update patient root
possible_data_dirs=dir(patient_root);
for b=1:length(possible_data_dirs)
    find_rec_dirs = regexp(possible_data_dirs(b).name,'adm_*','once');
     if ~isempty(find_rec_dirs)        
         patient_root = [possible_data_dirs(b).folder , '/', possible_data_dirs(b).name, '/'];
     end
end

root_files=dir(patient_root);
data_dirs=[];
dir_ct=0;

% Find data dirs
for a=1:length(root_files)
    dir_name = [root_files(a).folder , '/', root_files(a).name];

    find_rec_dirs = regexp(dir_name,'rec_*','once');
    if ~isempty(find_rec_dirs)
        dir_ct=dir_ct+1;
        data_dirs{dir_ct}=[dir_name, '/'];
    end
end


%% Get headers for patient
ieeg_fnames=[];
ieeg_dir=[];
ieeg_onset=[];
ieeg_offset=[];
ieeg_ct=0;
for a=1:length(data_dirs)
    fprintf('Searcing directory %s\n',data_dirs{a});
    
    %temp_files=dir(fullfile(data_dirs{a}));
    temp_files=dir(data_dirs{a});
    for b=1:length(temp_files)
        %         disp(temp_files(b).name);
        [stem, ext]=strtok(temp_files(b).name,'.');
        if strcmpi(ext,'.head')
            ieeg_ct=ieeg_ct+1;
            ieeg_fnames{ieeg_ct}=stem;
            ieeg_dir{ieeg_ct}=data_dirs{a};
            fname = fullfile(data_dirs{a},[ieeg_fnames{ieeg_ct} '.data']);
            temp_ieeg=bin_file(fname);
            ieeg_onset{ieeg_ct}=str2datetime(temp_ieeg.a_start_ts);
            %ieeg_onset{ieeg_ct}=datetime(temp_ieeg.a_start_ts,'InputFormat','yyyy-MM-dd hh:mm:ss');
            ieeg_offset{ieeg_ct}=str2datetime(temp_ieeg.a_stop_ts);
            %ieeg_onset{ieeg_ct}=temp_ieeg.a_start_ts;
            %ieeg_offset{ieeg_ct}=temp_ieeg.a_stop_ts;
        end
    end
end
fprintf('%d files for patient %s\n',ieeg_ct,patient_root);


%% Read list of seizures
in_fname=fullfile(patient_root,'seizurelist.txt');
csv_matrix=csv2Cell_EU(in_fname,9,4);

% Get list of files
szr_onsets=[];
szr_offsets=[];
szr_onsets_tpt=[];
szr_offsets_tpt=[];
szr_onsets_dt=[];
szr_offsets_dt=[];
n_szr=0;
for a=1:size(csv_matrix,1),
    if ~isempty(csv_matrix{a,1}),
        n_szr=n_szr+1;
        szr_onsets{n_szr}=csv_matrix{a,1}; % This is date-time of onset
        szr_onsets_tpt(n_szr)=str2num(csv_matrix{a,3}); % This must be the time point of onset in the data file
        %szr_onsets_dt(n_szr)=datetime(csv_matrix{a,1},'InputFormat','yyyy-MM-dd HH:mm:ss.s');
        szr_onsets_dt{n_szr}=str2datetime(csv_matrix{a,1});
        
        szr_offsets{n_szr}=csv_matrix{a,2}; % This is date-time of offset
        szr_offsets_tpt(n_szr)=str2num(csv_matrix{a,4}); % This must be the time point of offset in the data file
        %szr_offsets_dt(n_szr)=datetime(csv_matrix{a,2},'InputFormat','yyyy-MM-dd HH:mm:ss.s');
        szr_offsets_dt{n_szr}=str2datetime(csv_matrix{a,2});
    end
end


%% Loop over files and figure out which time points are ictal or not
for a=1:ieeg_ct,
    for b=1:n_szr,
        has_ictal=0;
        if  (ieeg_onset{a}<=szr_onsets_dt{b} && ieeg_offset{a}>=szr_onsets_dt{b}),
            % Szr begins in the file
            has_ictal=1;                
            temp_onset_sec=seconds(szr_onsets_dt{b}-ieeg_onset{a});
            if (ieeg_onset{a}<=szr_offsets_dt{b} && ieeg_offset{a}>=szr_offsets_dt{b}),
                % Szr also ends in the file
                temp_offset_sec=seconds(szr_offsets_dt{b}-ieeg_onset{a});
            else
                % Szr continues to the end of the file
                temp_offset_sec=NaN;
            end
        elseif ieeg_onset{a}<=szr_offsets_dt{b} && ieeg_offset{a}>=szr_offsets_dt{b},
            % Szr ends in file but began in previous file
            has_ictal=1;
            temp_onset_sec=0;
            temp_offset_sec=seconds(szr_offsets_dt{b}-ieeg_onset{a});
        end
        
        if has_ictal,
            fprintf('File %d contains szr %d\n',a,b);
            % Load header for the file
            temp_ieeg=bin_file(fullfile(ieeg_dir{a},[ieeg_fnames{a} '.data']));
            
            % Create a binary class vector
            labels=zeros(temp_ieeg.a_n_samples,1,'int8');
            temp_onset_tpt=1+round(temp_onset_sec*temp_ieeg.a_samp_freq);
            if isnan(temp_offset_sec),
                temp_offset_tpt=temp_ieeg.a_n_samples;
            else
                temp_offset_tpt=1+round(temp_offset_sec*temp_ieeg.a_samp_freq);
            end
            labels(temp_onset_tpt:temp_offset_tpt)=1;
            % Set according to onset/offset time
            fprintf('Start-Stop secs are: %d-%d\n',temp_onset_sec,temp_offset_sec);
            fprintf('Start-Stop tpts are: %d-%d\n',temp_onset_tpt,temp_offset_tpt);
            
            if abs(temp_onset_tpt-szr_onsets_tpt(b))>1,
                warning(sprintf('Onset of Szr#%d differs by more than one time point from seizurelist.txt',b));
            end
            if abs(temp_offset_tpt-szr_offsets_tpt(b))>1,
                warning(sprintf('Offset of Szr#%d differs by more than one time point from seizurelist.txt',b));
            end
            
            % Save labels
            out_fname=[ieeg_fnames{a} '.ictal'];
            fprintf('Saving ictal class vector to %s\n',out_fname);
            save(out_fname,'labels','temp_ieeg');
        end
    end
end



