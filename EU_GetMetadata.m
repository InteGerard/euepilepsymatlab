function [record_fnames, seizure_fnames, onset_samples, offset_samples, sc_seizure_fnames, sc_onset_samples, sc_offset_samples, ieeg_onset, ieeg_offset]= EU_GetMetadata(EU_data_path, patient_name)

%     global EU_data_path = '/media/gerard/166844A6160F0301/EU_DATA/inv/';

    %patient_name = 'pat_FR_1096';
    
    
    [current_dir,name,ext] = fileparts(mfilename('fullpath'))
    metadata_dir = [current_dir, '/metadata/'];
    
    subclin_fname = [metadata_dir, patient_name, '_subclinical_szrs.csv'];
    clin_fname = [metadata_dir, patient_name, '_clinical_szrs.csv'];
    

    subclin_foutname = [metadata_dir, patient_name, '_subclinical_szrs_detailed.csv'];
    clin_foutname = [metadata_dir, patient_name, '_clinical_szrs_detailed.csv'];
    
    
    patient_root=[EU_data_path, patient_name];

    % Find adm_ directory and update patient root
    possible_data_dirs=dir(patient_root);
    for b=1:length(possible_data_dirs)
        find_rec_dirs = regexp(possible_data_dirs(b).name,'adm_*','once');
         if ~isempty(find_rec_dirs)        
             patient_root = [patient_root, '/', possible_data_dirs(b).name, '/'];
         end
    end

    root_files=dir(patient_root);
    data_dirs=[];
    dir_ct=0;

    % Find data dirs
    for a=1:length(root_files)
        dir_name = [patient_root , '/', root_files(a).name];

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
                record_fnames{ieeg_ct}=fname;
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

    %% Read list of subclinical seizures
    in_fname=subclin_fname;
    csv_matrix=csv2Cell_EU(in_fname, ',');

    % Get list of files
    sc_szr_onsets=[];
    sc_szr_offsets=[];
    sc_szr_onsets_dt=[];
    sc_szr_offsets_dt=[];
    n_scszr=0;
    for a=1:size(csv_matrix,1),
        if ~isempty(csv_matrix{a,1}),
            n_scszr=n_scszr+1;
            sc_szr_onsets{n_scszr}=csv_matrix{a,1}; % This is date-time of onset
            %szr_onsets_dt(n_szr)=datetime(csv_matrix{a,1},'InputFormat','yyyy-MM-dd HH:mm:ss.s');
            sc_szr_onsets_dt{n_scszr}=str2datetime(csv_matrix{a,1});

            sc_szr_offsets{n_scszr}=csv_matrix{a,2}; % This is date-time of offset
            %szr_offsets_dt(n_szr)=datetime(csv_matrix{a,2},'InputFormat','yyyy-MM-dd HH:mm:ss.s');
            sc_szr_offsets_dt{n_scszr}=str2datetime(csv_matrix{a,2});
        end
    end  

    
    %% Loop over files and figure out which time points are subclinical or not
    incr = 0;
    for a=1:ieeg_ct,
        for b=1:n_scszr,
            has_ictal=0;
            if  (ieeg_onset{a}<=sc_szr_onsets_dt{b} && ieeg_offset{a}>=sc_szr_onsets_dt{b}),
                % Szr begins in the file
                has_ictal=1;                
                temp_onset_sec=seconds(sc_szr_onsets_dt{b}-ieeg_onset{a});
                if (ieeg_onset{a}<=sc_szr_offsets_dt{b} && ieeg_offset{a}>=sc_szr_offsets_dt{b}),
                    % Szr also ends in the file
                    temp_offset_sec=seconds(sc_szr_offsets_dt{b}-ieeg_onset{a});
                else
                    % Szr continues to the end of the file
                    temp_offset_sec=NaN;
                end
            elseif ieeg_onset{a}<=sc_szr_offsets_dt{b} && ieeg_offset{a}>=sc_szr_offsets_dt{b},
                % Szr ends in file but began in previous file
                has_ictal=1;
                temp_onset_sec=0;
                temp_offset_sec=seconds(sc_szr_offsets_dt{b}-ieeg_onset{a});
            end

            if has_ictal,
                %fprintf('File %d contains szr %d\n',a,b);

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
                %fprintf('Start-Stop secs are: %d-%d\n',temp_onset_sec,temp_offset_sec);
                %fprintf('Start-Stop tpts are: %d-%d\n',temp_onset_tpt,temp_offset_tpt);

                incr = incr+1;
                
                sc_seizure_fidx{incr} = a;
                sc_seizure_fnames{incr} = record_fnames{a};
                sc_onset_dts{incr} = datestr(sc_szr_onsets_dt{b});
                sc_offset_dts{incr} = datestr(sc_szr_offsets_dt{b});
                sc_onset_samples{incr} =  temp_onset_tpt;
                sc_offset_samples{incr} = temp_offset_tpt;

                % Save labels
%                 out_fname=[ieeg_fnames{a} '.ictal'];
%                 fprintf('Saving ictal class vector to %s\n',out_fname);
%                 save(out_fname,'labels','temp_ieeg');
            end
        end
    end  
    
    scsz_data = [sc_seizure_fidx ; sc_seizure_fnames; sc_onset_dts; sc_offset_dts; sc_onset_samples; sc_offset_samples]';
    
     
    
    
    %% Read list of seizures
    in_fname=clin_fname;
    csv_matrix=csv2Cell_EU(in_fname, ',');

    % Get list of files
    szr_onsets=[];
    szr_offsets=[];
    %szr_onsets_tpt=[];
    %szr_offsets_tpt=[];
    szr_onsets_dt=[];
    szr_offsets_dt=[];
    n_szr=0;
    for a=1:size(csv_matrix,1),
        if ~isempty(csv_matrix{a,1}),
            n_szr=n_szr+1;
            szr_onsets{n_szr}=csv_matrix{a,1}; % This is date-time of onset
            %szr_onsets_tpt(n_szr)=str2num(csv_matrix{a,3}); % This must be the time point of onset in the data file
            %szr_onsets_dt(n_szr)=datetime(csv_matrix{a,1},'InputFormat','yyyy-MM-dd HH:mm:ss.s');
            szr_onsets_dt{n_szr}=str2datetime(csv_matrix{a,1});

            szr_offsets{n_szr}=csv_matrix{a,2}; % This is date-time of offset
            %szr_offsets_tpt(n_szr)=str2num(csv_matrix{a,4}); % This must be the time point of offset in the data file
            %szr_offsets_dt(n_szr)=datetime(csv_matrix{a,2},'InputFormat','yyyy-MM-dd HH:mm:ss.s');
            szr_offsets_dt{n_szr}=str2datetime(csv_matrix{a,2});
        end
    end

    %% Loop over files and figure out which time points are ictal or not
    incr = 0;
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
                %fprintf('File %d contains szr %d\n',a,b);

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
                %fprintf('Start-Stop secs are: %d-%d\n',temp_onset_sec,temp_offset_sec);
                %fprintf('Start-Stop tpts are: %d-%d\n',temp_onset_tpt,temp_offset_tpt);

%                 if abs(temp_onset_tpt-szr_onsets_tpt(b))>1,
%                     warning(sprintf('Onset of Szr#%d differs by more than one time point from seizurelist.txt',b));
%                 end
%                 if abs(temp_offset_tpt-szr_offsets_tpt(b))>1,
%                     warning(sprintf('Offset of Szr#%d differs by more than one time point from seizurelist.txt',b));
%                 end

                incr = incr+1;
                
                seizure_fidx{incr} = a;
                seizure_fnames{incr} = record_fnames{a};
                onset_dts{incr} = datestr(szr_onsets_dt{b});
                offset_dts{incr} = datestr(szr_offsets_dt{b});
                onset_samples{incr} =  temp_onset_tpt;
                offset_samples{incr} = temp_offset_tpt;

                % Save labels
%                 out_fname=[ieeg_fnames{a} '.ictal'];
%                 fprintf('Saving ictal class vector to %s\n',out_fname);
%                 save(out_fname,'labels','temp_ieeg');
            end
        end
    end

    

    sz_data = [seizure_fidx ; seizure_fnames; onset_dts; offset_dts; onset_samples; offset_samples]';


    % Write out detailed files for use in python later
    cell2csv(clin_foutname, sz_data, ',');
    cell2csv(subclin_foutname, scsz_data, ',');
    
    
end