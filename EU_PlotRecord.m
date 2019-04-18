% Read data from the EU database files. Use metadata from http://epilepsiae.uniklinik-freiburg.de
% to annotate these files with the clinical and subclincal seizure onset and end times.
% GOL 2017 - University of Toronto
function [ ] = EU_PlotRecord(record, image_dir)

    signal = record.signal;
    Fs = record.Fs;
    channels = record.channels;
    selected_record = record.selected_record;
    patient_id = record.patient_id;
    sample_sz_onset = record.sample_sz_onset;
    sample_sz_end = record.sample_sz_end;
    sample_scsz_onset = record.sample_scsz_onset;
    sample_scsz_end = record.sample_scsz_end;
    sz_idx = record.isSZ;
    scsz_idx = record.isSCSZ;
    ieeg_onset = record.ieeg_onset;
    ieeg_offset = record.ieeg_offset;
    record_idx = record.record_idx;

    %% Plot Ictal Activity  
    fig = figure(2);
    set(gcf,'Position', get(0,'Screensize')); % Maximize figure.
    set(gca,'LooseInset',get(gca,'TightInset'));
    %set(gcf, 'Visible', 'off')
    clf(fig)

    norm_factor = max(max(abs(signal))); % max absolute value in all time series
    n_channels = size(signal,2);
    for plot_channel_num = 1:n_channels

        x = signal(:,plot_channel_num);
        x = x - mean(x);
        channel_sig = x / norm_factor;

        channel_sig = channel_sig + (n_channels - plot_channel_num + 1); % adjust plot
        %channel_sig = 2*(x - min(x))/(max(x) - min(x)) - 1; % Normalize

        siglen = size(signal,1);
        tscale = linspace(0,siglen./Fs,siglen)';
        plot(tscale,channel_sig, 'k');
        hold on;
        axis([0 max(tscale) 0 n_channels+1]);

        if (sz_idx) 
            for i = 1:length(sample_sz_onset)
                sample_onset = sample_sz_onset{i};
                sample_end = sample_sz_end{i};

                tscale = linspace(sample_onset/Fs,sample_end./Fs,sample_end-sample_onset)';
                seizure_signal = channel_sig(sample_onset:sample_end-1,:);

                plt = plot(tscale,seizure_signal,'r');
            end
        end

        if (scsz_idx) 
            for i = 1:length(sample_scsz_onset)
                sample_onset = sample_scsz_onset{i};
                sample_end = sample_scsz_end{i};

                tscale = linspace(sample_onset/Fs,sample_end./Fs,sample_end-sample_onset)';
                seizure_signal = channel_sig(sample_onset:sample_end-1,:);

                plot(tscale,seizure_signal,'g');
            end                    
        end
    end
    [a,b] = regexp(patient_id,['\d+']);
    pat_num = patient_id(a:b);
    title(sprintf( 'Patient %s, Record %d: %s', pat_num, record_idx, datestr(ieeg_onset) ) );
    xlabel('Time (Seconds)');
    ylabel('Channel Name');
    set(gca,'YTick', 1:1:n_channels, 'YTickLabels', channels);
    saveas(fig,sprintf( '%s/raw_%s_%i.bmp', image_dir, patient_id, record_idx ) );
end       