%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Framework for feature extraction with EU epilepsy database
% Gerard O'Leary - University of Toronto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ done ] = EU_GetPatientRecord (patient_id, EU_data_path, output_dir)

    done = 0;

    % Wait for x seconds before executing
    %pause(7200);
    % profile -memory on;
    % setpref('profiler','showJitLines',1);
    % profile viewer

    addpath(genpath('../../NURIP/Feature_Extraction/online'));
    addpath(genpath('../../EU'));

    %EU_data_path = '/media/gerard/My Book/EU/inv/';
    %EU_data_path = '/media/gerard/166844A6160F0301/EU_DATA/inv/';

    %% Specify patient of records to process
    % patient_id = 'pat_FR_1096' 'pat_FR_970' 'pat_FR_548' 'pat_FR_1125' 'pat_FR_442';
    %patient_id = 'pat_FR_1096';




    %% Set up configuration and fiele access paths
    [channels, se_alphas, se_bands, plv_alphas, plv_bands, downsample_rate, records_to_process] = PRJ_LoadConfig(patient_id);

    %% DEGUG
    %records_to_process = 14;



    %output_dir = './Feature_Extraction/EU_OUTPUT/';
    %output_dir = '/media/gerard/SSD_256GB/EU_Processing/EDMSE_PLV/';
    mkdir(output_dir);

    output_dir = [output_dir, patient_id];
    mkdir(output_dir);

    image_dir = [output_dir, '/images/'];
    mkdir(image_dir);

    %% Process patient records
    dataset_ds_rate = 10; % downsample rate for the final extracted features - 4 should be ok for the EDM-SE feature as it's recursive.

    for record_idx = records_to_process


        [record] = EU_GetRecord(EU_data_path, patient_id, record_idx, channels, downsample_rate);
        EU_PlotRecord(record, image_dir);

        %signal = signal(1:1000,:);


        %img_fname = sprintf( '%s/spec_%s_%i.bmp', spec_image_dir, patient_id, record_idx );
        %ENERGY_PlotSpectrogram ( signal, Fs, sample_sz_onset, img_fname );

        tstart = datestr(now);

        %signal = signal(700000:800000,:);
        [se_concat, se_means, se_sds] = run_EDMSE_online( record.signal, se_bands, se_alphas, channels);
        %[plv_concat, plv_means, plv_sds] = PLV_ExtractAllCombs( signal, plv_bands, plv_alphas, channels);

        tend = datestr(now);

        fprintf('\nFinished Extracting index: %i \n', record_idx);
        fprintf('Start: %s\n',tstart);
        fprintf('End: %s\n',tend);
        

        %% Downsample before saving
        se_concat = downsample(se_concat,dataset_ds_rate);
%         plv_concat = downsample(plv_concat,dataset_ds_rate);
        for i = 1:length(record.sample_sz_onset)   
            record.sample_sz_onset{i} = floor(record.sample_sz_onset{i}/dataset_ds_rate);
            record.sample_sz_end{i} = floor(record.sample_sz_end{i}/dataset_ds_rate);
        end

        for i = 1:length(sample_scsz_onset)   
            record.sample_scsz_onset{i} = floor(sample_scsz_onset{i}/dataset_ds_rate);
            record.sample_scsz_end{i} = floor(sample_scsz_end{i}/dataset_ds_rate);
        end
                

        %% Save and plot
        fname = [ output_dir, '/EDMSE_', patient_id , '_', num2str(record_idx, '%03d'), '.mat'];
        
        signal = record.signal;
        Fs = record.Fs;
        record_fname = record.record_fname;
        sample_sz_onset = record.sample_sz_onset;
        sample_sz_end = record.sample_sz_end;
        sample_scsz_onset = record.sample_scsz_onset ;
        sample_scsz_end = record.sample_scsz_end;
    
        save(fname, ...
            'patient_id', ...
            'record_idx', ...
            'record_fname', ...
            'Fs', ...
            'sample_sz_onset', ...
            'sample_sz_end', ...
            'sample_scsz_onset', ...
            'sample_scsz_end', ...
            'signal',...
            'se_means', ...
            'se_sds', ...        
            'se_concat', ...
            'se_bands', ...
            'se_alphas', ...        
            'dataset_ds_rate', ...
            'tstart', ...
            'tend', ...
            '-v7.3' ); % Use when filename is stored in a variable
            %commented to avoid compression % uncommented 7.3 as it's
            %necessary to use the hdf5 library in Python.

        fig = figure(10);
        clf(fig);
        features = [se_concat];
        features = zscore(features);
        imagesc(features');
        caxis([-3, 3]);
        title('Feature Space');
        colorbar;
        img_fname = sprintf( '%s/feature_space_%s_%i.bmp', image_dir, patient_id, record_idx);
        xlabel('Time (samples)')
        ylabel('EDMSE');
        saveas(fig,img_fname);

                %se_concat = zscore(se_concat,1,1);


        %img_fname = sprintf( '%s/edmse_%s_%i.bmp', edmse_image_dir, patient_id, record_idx );

    %     if(EN_PLOT)
    %         PlotEDMSE ( se_concat, Fs, sample_sz_onset, sample_sz_end, img_fname );
    %     end

    %     min = -1
    %     max = 1
    %     t1 = se_concat(se_concat>max)=max;
    %     t2 = se_concat(se_concat<min)=min;
        
        

    end
    
    done = 1;
end