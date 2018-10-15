%% FlexVal eventMark analysis -- IR trigger + Extender via EPOC
%{
Author: Nik Willimas
Date: August 2018
Purpose: Process, plot, and provide descriptive statistics for
event-marking hardware

Glossary:
    data.johnson_eeg --> array of the eeg data from the johnson infra-red trigger
    data.johnson_diff --> array of the derivative (change) from sample to sample
        in the johnson eeg waveform
    data.johnson_channel --> array of derived johnson events
    data.extender_channel --> array of extender events

    NOTE: The following time variables are represented in both samples
        and ms (with _ms extenstion).
    times.extender_events --> array of the actual times at which parallel port
        events occured.
    times.johnson_events --> array of the actual times at which the derived
        johnson events are calculated to have occured.
    times.johnson_latency --> array of the delay of each derived johnson
        event relative to each parallel port events
    times.johnson_interval --> array of the time between each johnson event
    times.pport_interval --> array of the time between each parallel port event
%}

%% Parameters and Inputs
%Toggles
r.import = 1;
r.filter = 1;
r.plot = 1;
r.write = 0;
r.data_set = 5; %which data set to use
%   1 = 1000 matTones --> extender + EPOC 256hz
%   2 = 1000 fireface clicks --> extender + EPOC 256hz
%   3 = 1000 matTones --> extender + Flex 128hz
%   4 = 1000 fireface clicks --> extender + Flex 128hz
%   5 = 1000 matTones --> extender + EPOC 128hz
%   6 = 1000 fireface clicks --> extender + EPOC 128hz
%   7 = 1000 fireface clicks --> SD mode extender + EPOC 128hz
%   8 = 1000 psychopyAPI tones --> API + EPOC 256hz

p.epoc_samplingRate = 253.7;
p.epoc_samplingRate_128 = 126.9;
p.flex_samplingRate = 129.1;
p.extenderSD_samplingRate = 126.9;

%Switch for different parameters depending on data set
switch r.data_set 
    case 1
        %%Parameters
        p.name = 'EPOC256_IR_Extender_matTones.csv'; %name for export files
        p.jDiff_multiplier = 10; %magnify johnson derivative
        p.extenderEvent_size = 200; %size to make parallel port event tics for graphing
        p.johnson_threshold = 20; %set johnson derivative threshold to identify johnson events
        p.jEventMark_size = 200; %how big to make the event mark (for graphing)
        p.jOnset_multiplier = 50; %multiply the onset value by something for graphing
        p.x_range = 14200:14600; % what to show on the x-axis for graphs
        p.plot_title = '1000 Matlab Tones - Extender w/ EPOC (256hz)';
        p.sampling_rate = p.epoc_samplingRate;
    case 2
        %%Parameters
        p.name = 'EPOC256_IR_Extender_fireface.csv'; %name for export files
        p.jDiff_multiplier = 7; %magnify johnson derivative
        p.extenderEvent_size = 200; %size to make parallel port event tics for graphing
        p.johnson_threshold = 28; %set johnson derivative threshold to identify johnson events
        p.jEventMark_size = 200; %how big to make the event mark (for graphing)
        p.jOnset_multiplier = 50; %multiply the onset value by something for graphing
        p.x_range = 14250:14700; % what to show on the x-axis for graphs
        p.plot_title = '1000 Fireface Clicks - Extender w/ EPOC (256hz)';
        p.sampling_rate = p.epoc_samplingRate;
    case 3
        %%Parameters
        p.name = 'Flex128_IR_Extender_matTones.csv'; %name for export files
        p.jDiff_multiplier = 20; %magnify johnson derivative
        p.extenderEvent_size = 100; %size to make parallel port event tics for graphing
        p.johnson_threshold = 22; %set johnson derivative threshold to identify johnson events
        p.jEventMark_size = 100; %how big to make the event mark (for graphing)
        p.jOnset_multiplier = 50; %multiply the onset value by something for graphing
        p.x_range = 14250:14500; % what to show on the x-axis for graphs
        p.plot_title = '1000 Matlab Tones - Extender w/ Flex';
        p.sampling_rate = p.flex_samplingRate;
    case 4
        %%Parameters
        p.name = 'Flex128_IR_Extender_fireface.csv'; %name for export files
        p.jDiff_multiplier = 10; %magnify johnson derivative
        p.extenderEvent_size = 10; %size to make parallel port event tics for graphing
        p.johnson_threshold = 22; %set johnson derivative threshold to identify johnson events
        p.jEventMark_size = 10; %how big to make the event mark (for graphing)
        p.jOnset_multiplier = 10; %multiply the onset value by something for graphing
        p.x_range = 14200:14450; % what to show on the x-axis for graphs
        p.plot_title = '1000 Fireface Clicks - Extender w/ Flex';
        p.sampling_rate = p.flex_samplingRate;
    case 5
        %%Parameters
        p.name = 'EPOC128_IR_Extender_matTones.csv'; %name for export files
        p.jDiff_multiplier = 7; %magnify johnson derivative
        p.extenderEvent_size = 200; %size to make parallel port event tics for graphing
        p.johnson_threshold = 52; %set johnson derivative threshold to identify johnson events
        p.jEventMark_size = 200; %how big to make the event mark (for graphing)
        p.jOnset_multiplier = 57; %multiply the onset value by something for graphing
        p.x_range = 14200:14600; % what to show on the x-axis for graphs
        p.plot_title = '1000 Matlab Tones - Extender w/ EPOC (128hz)';
        p.sampling_rate = p.epoc_samplingRate_128;
    case 6
        %%Parameters
        p.name = 'EPOC128_IR_Extender_fireface.csv'; %name for export files
        p.jDiff_multiplier = 7; %magnify johnson derivative
        p.extenderEvent_size = 200; %size to make parallel port event tics for graphing
        p.johnson_threshold = 51; %set johnson derivative threshold to identify johnson events
        p.jEventMark_size = 200; %how big to make the event mark (for graphing)
        p.jOnset_multiplier = 50; %multiply the onset value by something for graphing
        p.x_range = 14250:14700; % what to show on the x-axis for graphs
        p.plot_title = '1000 Fireface Clicks - Extender w/ EPOC (128hz)';
        p.sampling_rate = p.epoc_samplingRate_128;
    case 7
        %%Parameters
        p.name = 'EPOC128_ExtenderSD_fireface';
        p.extenderEvent_size = 200; %size to make parallel port event tics for graphing
        p.x_range = 14200:14600; % what to show on the x-axis for graphs
        p.plot_title = '1000 Fireface Clicks - Extender SD mode (128hz)';
        p.sampling_rate = p.epoc_samplingRate_128;
    case 8
        %%Parameters
        p.name = 'EPOC256_psychopyAPI'; %name for export files
        p.jDiff_multiplier = 7; %magnify johnson derivative
        p.event_size = 200; %size to make parallel port event tics for graphing
        p.jOnset_multiplier = 50; %multiply the onset value by something for graphing
        p.x_range = 14250:14700; % what to show on the x-axis for graphs
        p.plot_title = '1000 psychoppyAPI tones - API w/ EPOC (256hz)';
        p.sampling_rate = p.epoc_samplingRate;
end

%% Data export
p.dir_out = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/matlab_ouput';
p.fullfile_out = fullfile(p.dir_out,p.name);

%% Data Import
if r.import
    %file locations
    switch r.data_set
        case 1
            %1000 trials - matTones - Extender w/ EPOC
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/extender_vs_IR/epoc_256hz/matTones/FlexVal_IR_Extender_EPOC_2018.08.13_08.41.48.edf';        
        case 2
            %1000 trials - fireface clicks - Extender w/ EPOC
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/extender_vs_IR/epoc_256hz/fireface/FlexVal_eventMark_Fireface_IR_Extender_EPOC_1000_2018.08.15_12.41.41.edf';
        case 3
            %1000 trials - matTones - Extender w/ Flex
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/extender_vs_IR/flex/matTones/FlexVal_IR_Extender_Flex_MatlabTones1000_2018.08.29_10.25.40.edf';
        case 4
            %1000 trials - fireface clicks - Extender w/ Flex
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/extender_vs_IR/flex/fireface/FlexVal_eventMark_IR_Extender_Fireface_Flex_1000_2018.08.29_09.52.06.edf';
        case 5
            %1000 trials - matTones - Extender w/ EPOC
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/extender_vs_IR/epoc_128hz/matTones/FlexVal_eventMark_IR_Extender_EPOC_128hz_matTones_2018.09.05_08.29.04.edf';
        case 6
            %1000 trials - fireface clicks - Extender w/ EPOC
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/extender_vs_IR/epoc_128hz/fireface/FlexVal_eventMark_IR_Extender_EPOC_128hz_fireface1000_2018.09.05_10.38.57.edf';
        case 7
            %1000 trials - fireface clicks - SD Extender w/ EPOC 128
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/extender_SDmode/2018_09_05_11_07_59__001_eeg.csv';
        case 8
            %1000 trials - psychopy tones - API w/ EPOC 256
            currentData = '/Users/mq43504108/Macquarie University/Nicholas Badcock - NikWilliams/Projects/eventMark/data/psychopyAPI/eventMark_psychopyAPI_IR_1000_2018.10.08_09.19.30.edf';
    end
    
    %bring in the data
    
    
    switch r.data_set
        case {1,2,3,4,5,6,8}
            EEG = pop_biosig(currentData); %can change path
            %rename channels into something sensible
            emotiv.channels = [];
            for i = 1 : numel(EEG.chanlocs); emotiv.channels{i} = EEG.chanlocs(i).labels; end
            
            %get EPOC data stream channels of interest
            switch r.data_set
                case {1,2,5,6}
                    emotiv.chans_wanted = {'T7', 'T8', 'MARKER_HARDWARE'};
                    emotiv.chans_wanted_index = ismember(emotiv.channels, emotiv.chans_wanted);
                case {3,4}
                    emotiv.chans_wanted = {'TP7', 'TP8', 'MARKER_HARDWARE'};
                    emotiv.chans_wanted_index = ismember(emotiv.channels, emotiv.chans_wanted);
                case 8
                    emotiv.chans_wanted = {'T7', 'T8', 'MARKER_HARDWARE'};
                    emotiv.chans_wanted_index = ismember(emotiv.channels, emotiv.chans_wanted);
            end
            
            %% Processing
            %get all the data in one matrix
            data.matrix = EEG.data;
            %just get the stuff we're interested in
            data.matrix_small = data.matrix(emotiv.chans_wanted_index,:);
            %Johnson IR trigger data w/ geometric transformation so that the y-axis is sensible
            data.johnson_eeg = data.matrix_small(1,:) - 4200;
            %Extender events with a multiplier so that they show up on plots
            data.extender_channel = data.matrix_small(3,:)* p.extenderEvent_size;
            
            if r.filter
                %filter the Johnson signal
                low_pass=40.0;
                high_pass=0.5;
                fs=p.sampling_rate;
                [b, a] = butter(2, [high_pass, low_pass] / fs);
                data.johnson_eeg = filtfilt(b,a,double(data.johnson_eeg));
            end
            
            %{
calculate the derivative of the EEG waveform. Essentially this is n -
(n+1) at each sample, n. This will give an array where large values indicate a sudden change in
voltage
            %}
            %data.johnson_diff = [0,diff(data.johnson_eeg)]* p.jDiff_multiplier;
            data.johnson_diff = [0,diff(data.johnson_eeg)];
            
            %calculate the onset of the signal
            %data.johnson_onset = (data.johnson_diff > p.johnson_threshold) * p.jOnset_multiplier;
            data.johnson_onset = (data.johnson_diff > p.johnson_threshold);
            
            %{
For each epoch, take the first instance of a johnson_onset, set it equal
to something non-zero, and set all the other values equal to zero. This is
done because there are multiple places within each epoch that the
johnson_diff value exceeds the threshold (we only care about the first).
            %}
            
            %how big of a slice to take around each pport event (in ms)
            extender_epoc_interval = [-100 100];
            
            times.extender_events = find(data.extender_channel);
            %create an array which contains the samples at which extender marked an
            %event (we use this to index the eeg data)
            times.extender_events = times.extender_events(2:length(times.extender_events)); %chuck out the first event (it's when recording starts)
            
            %make a new channel which contains the derived johnson events
            data.johnson_channel = zeros(1, length(data.johnson_onset));
            
            for i = 1 : numel(times.extender_events)
                current_extender_event = times.extender_events(i);
                %current_interval_end = current_extender_event + extender_epoc_interval;
                %current_interval = [current_extender_event current_interval_end];
                current_interval = current_extender_event + extender_epoc_interval;
                current_johnson_onsets = data.johnson_onset(current_interval(1):current_interval(2));
                current_johnson_first = find(current_johnson_onsets ~=0, 1, 'first');
                new_johnson_onsets = zeros(1, length(current_johnson_onsets));
                new_johnson_onsets(current_johnson_first) = p.jEventMark_size;
                data.johnson_channel(current_interval(1):current_interval(2)) = new_johnson_onsets;
            end
            
            %%
            %Timing
            %at which samples are the johnson events?
            times.johnson_events = find(data.johnson_channel);
            % how much slower is the johnson trigger relative to the parallel port
            % event?
            
            %Make a nice round 1000 events
            times.extender_events = times.extender_events(1:1000);
            times.johnson_events = times.johnson_events(1:1000);
            
            %THESE ARE IN SAMPLES
            times.johnson_latency = times.johnson_events - times.extender_events;
            
            %what is the interval between johnson events?
            times.johnson_interval = diff(times.johnson_events);
            
            %what is the interval between parallel port events?
            times.extender_interval = diff(times.extender_events);
                      
            data.times = (1:size(data.matrix,2))/p.sampling_rate;
            
            times.extender_events_ms = (times.extender_events/p.sampling_rate)*1000;
            times.johnson_events_ms = (times.johnson_events/p.sampling_rate)*1000;
            times.johnson_latency_ms = (times.johnson_events_ms - times.extender_events_ms);
            times.johnson_interval_ms = diff(times.johnson_events_ms);
            times.extender_interval_ms = diff(times.extender_events_ms);
            
            %%Descriptives
            
            data.times_fields = {'johnson_latency', 'johnson_interval', 'extender_interval'};
            
            fprintf('\n%s\n%s\n', p.plot_title,'Timing results (in samples):')
            for i = 1:numel(data.times_fields)
                current_mean = mean(times.(data.times_fields{i}));
                current_std = std(times.(data.times_fields{i}));
                current_field = data.times_fields{i};
                fprintf('\n%s%s%s%.2f%s%.2f\n', current_field,':', ' M = ', current_mean, ', SD = ', current_std)
            end
            
            data.times_fields_ms = {'johnson_latency_ms', 'johnson_interval_ms', 'extender_interval_ms'};
            
            fprintf('\n%s\n%s\n', p.plot_title,'Timing results (in ms):')
            for i = 1:numel(data.times_fields)
                current_mean = mean(times.(data.times_fields_ms{i}));
                current_std = std(times.(data.times_fields_ms{i}));
                current_field = data.times_fields_ms{i};
                plot_string{i} = sprintf('%s%s%s%.2f%s%.2f', current_field,':', ' M = ', current_mean, ', SD = ', current_std);
                fprintf('\n%s%s%s%.2f%s%.2f\n', current_field,':', ' M = ', current_mean, ', SD = ', current_std)
            end
            %%Perform Bartletts Test and print to console
            data_agg_ms = [times.johnson_interval_ms(:), times.extender_interval_ms(:)];
            [b,stats] = vartestn(data_agg_ms);
            
            bartletts_string = sprintf('\n%s\n%s%.2f%s%.3f\n','Bartletts test of equality of variance:','Bartletts statistic = ',...
                stats.chisqstat, ', p = ', b);
            disp(bartletts_string);
            
            %% Array export (for plotting in R)
            if r.write
                array_new = [times.extender_interval_ms; times.johnson_interval_ms]';
                table_out = array2table(array_new, 'Variablenames',...
                    {'extender_interval','johnson_interval'});
                writetable(table_out,p.fullfile_out);
            end
            
            %% Plotting
            if r.plot
                close; figure; hold;
                data.plot_fields = {'johnson_eeg','johnson_diff','johnson_channel','extender_channel',};
                data.legend_names = {};
                data.colours = {'b','g','r','k'};
                for i = 1 : numel(data.plot_fields)
                    plot(data.times(p.x_range),data.(data.plot_fields{i})(p.x_range),data.colours{i});
                    data.legend_names{i} = strrep(data.plot_fields{i},'_',' ');
                end
                title(p.plot_title);
                legend(data.legend_names);
                annotation('textbox', [.5 .01 .3 .3],'String',plot_string,'FitBoxToText','on', 'Interpreter','none');
                
                interval_hist = figure;
                
                %calculate the largest deviation from 1000 in the interval arrays so we can have the same bins for
                %both histograms
                hist_edges = [1000 - [min([times.johnson_interval_ms(:);times.extender_interval_ms(:)])]...
                    [max([times.johnson_interval_ms(:);times.extender_interval_ms(:)])] - 1000];
                hist_edges_max = ceil(max(hist_edges));
                hist_bins = [1000 - hist_edges_max : 1000 + hist_edges_max;];
                
                %calculate the maximum frequency so y-axis can be set accordingly
                yMax = max([histcounts(times.johnson_interval_ms) histcounts(times.extender_interval_ms)]);
                
                %make the histograms
                subplot(1,2,1)
                johnson_hist = histogram(times.johnson_interval_ms, hist_bins);
                title('Johnson Event Intervals');
                xlabel('Milliseconds');
                ylim([0 yMax]);
                t1 = annotation('textbox', [.2 .9 0 0],'String',plot_string{2},'FitBoxToText','on',...
                    'Interpreter','none','EdgeColor','none');
                t1.FontSize = 10;
                subplot(1,2,2)
                extender_hist = histogram(times.extender_interval_ms, hist_bins);
                title('Extender Event Intervals');
                xlabel('Milliseconds');
                ylim([0 yMax]);
                t2 = annotation('textbox', [.6 .9 0 0],'String',plot_string{3},'FitBoxToText','on',...
                    'Interpreter','none','EdgeColor','none');
                t2.FontSize = 10;
            end
        case 7
            eeg.table = readtable(currentData); %bring in data
            data.extender_channel = table2array(eeg.table(:,'MARKER_HARDWARE'));
            times.extender_events = find(data.extender_channel>0);
            times.extender_interval = diff(times.extender_events);
            times.extender_interval_ms = (times.extender_interval/p.extenderSD_samplingRate)*1000;
            missing_events_index = find(times.extender_interval_ms>1100); %see if there were any missed events
            missing_events_count = numel(missing_events_index);
            if numel(missing_events_count > 0)    %if there were, remove those intervals
                times.extender_interval(missing_events) = [];
                times.extender_interval_ms(missing_events) = [];
            end
            times.extender_interval = times.extender_interval(1:999);      %make it a nice round 1000 events
            times.extender_interval_ms = times.extender_interval_ms(1:1000);
            fprintf('\n%s\n%s%.2f\n%s%.2f\n','Extender interval (samples) ', 'M = ', mean(times.extender_interval),...
                'SD = ', std(times.extender_interval));
            fprintf('\n%s\n%s%.2f\n%s%.2f\n','Extender interval (ms) ', 'M = ', mean(times.extender_interval_ms),...
                'SD = ', std(times.extender_interval_ms));
            
    end
end

