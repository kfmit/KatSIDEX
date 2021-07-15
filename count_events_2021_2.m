% count events v. time, create detections w/ number of nodes.
close all
clear all
clc

%% A block for saving to a .mat, already done

%mat_out = '/media/efischell/Samsung_T5/events_metadata/';
mydir1='*';
nameout='pre0324_retest_pt2';
%prefix = ['/media/efischell/Samsung_T5/2021_pre0324/' mydir1 '/'];
%mat_out = '/media/efischell/Samsung_T5/events_metadata_new2/';
mat_out = '/home/kfung/Downloads/SIDEx/data21';
%prefix = ['/media/efischell/Samsung_T5/2021_pre0324/' mydir1 '/'];
prefix = ['/home/kfung/Downloads/SIDEx/data21/SIDEx_2021_field_data/' mydir1 '/']
min_select=9; % minimum number of detects needed to plot event; count event no matter what!

%mkdir(mat_out)

%% for all of the data in prefix/*/*.txt, detect events.
%directory = dir([prefix 'Sidex_20210324T1718*'])%%404*.txt']);
directory = dir([prefix 'Sidex_2021*.txt'])
FS=1000; % sample rate
ch_plot=1:32;%channels to plot
N=size(directory);

%% for each file, get timestamp:
epochs=[];
set(0,'DefaultFigureVisible','off')
for ii=1:N
    [epochs(ii),datetimes]=timeFromFilename(directory(ii).name);
end


e_start_ind=1;
e_start = epochs(e_start_ind);

data=[];
t=[];
event_times_all=[];
event_numdetect_all=[];
filenames_txt={};
start_iterate=1;
for ii=14788:N
    cur_time = epochs(ii);
    elapsed = cur_time-e_start;
    
%     filename = [directory(ii).folder '/' directory(ii).name]
    try 
        M = dlmread(filename, ',', 2, 0);
        cur_timevec=cur_time+(1:size(M,1))/FS;
    catch
        disp('file error')
        M=[];
        cur_timevec=[];
    end
    t=cur_timevec;
    if size(M,2)==32
        % then we have the data! Look for events.
        data=(M-mean(M)) / (2^24/2.5);%M;% do one minute at a time
        [b,a]=butter(6,200/(1000/2),'low');
        [c,d]=butter(6,1/(1000/2),'high');
        %plot(data)
        %pause
        %% first, find events
        event_ts=[];
        event_ts_all=[];
        % make directory for event .mat files:
                string_date = directory(ii).name(7:14);
                mat_loc = [mat_out 'events_' string_date '/'];
                mkdir(mat_loc)
                
        % use max channel to exclude cross-talk:
        ch_check = [1:15 19:31]; 
        
        [val,ch_corr]=max(max(data(:,ch_check)));
        data_filt_comp = filter(b,a,data(:,ch_check(ch_corr)));
        
        data_filt_comp = filter(c,d,data_filt_comp);
        background_filt =  movmedian(abs(data_filt_comp),FS*2);
        data_filt_comp_norm = (abs(data_filt_comp))./background_filt;
        
        % find the peaks:
        [~,pk_list_comp] = findpeaks(data_filt_comp_norm,'MinPeakProminence',10,'MinPeakDistance',1000);
        if length(pk_list_comp)==0
           pk_list_comp=[-1000]; 
        end
        for ch=ch_plot
            data_filt = data(:,ch);% filter(b,a,data(:,ch));
            %data_filt = filter(c,d,data_filt);
            
            % on filtered data, subtract out background noise;
            background = movmedian(abs(data_filt),FS*60);

            data_norm_log = (abs(data_filt))./background;
            
            if max(data_norm_log) > 3
                
            if ch ~= 16 && ch ~=32
                [~,pk_list] = findpeaks(data_norm_log,'MinPeakProminence',8,'MinPeakDistance',1000);
                findpeaks(data_norm_log,'MinPeakProminence',8,'MinPeakDistance',1000);
                new_pks=[];
                % for each event, check for correlation:
                for lll=1:length(pk_list)
                    start_val = max(1,pk_list(lll)-200);
                    end_val = min(length(data_norm_log),pk_list(lll)+200);
                [c,lags]=xcorr(data_norm_log(start_val:end_val), data_filt_comp_norm(start_val:end_val));
                [val,ind]=max(c);
                if min(abs(pk_list(lll)-pk_list_comp))>5 && abs(lags(ind))>5 && abs(data_filt_comp(pk_list(lll)))<0.1 && abs(data_filt(pk_list(lll)))>max(max(data(:,ch_check)))/100 % must be significant phase shift
                    lags(ind);
                    
                    new_pks=[new_pks; pk_list(lll)];
                    
                   % plot(data_norm_log(start_val:end_val));hold on;plot(data_filt_comp_norm(start_val:end_val))
                   % pause
                   % plot(lags,c)
                %pause
                end
                
                end
                event_ts=[event_ts; new_pks];
                event_ts_all = [event_ts_all; pk_list];
            end
            end
            
        end
        if length(event_ts)>0
            
        %% check for detections: first, get unique
        % we'll look at 2 s at a time, so any events w/in 2 s are assumed similar
        bin_size = 2; % bin size in seconds
        event_ts_round = ceil((event_ts/FS)/bin_size)*bin_size*FS;
        event_ts_all_round=ceil((event_ts_all/FS)/bin_size)*bin_size*FS
        [temp,inds_unique] = unique(event_ts_round)
        
        event_ts_unique = event_ts(inds_unique);
        
        %% for each event, save number of detects, if > min_select, plot
        for ll=1:length(event_ts_unique)
            
            %%
            time_idx = event_ts_unique(ll);
            time_idx_round = ceil((time_idx/FS)/bin_size)*bin_size*FS;
            num_detects = length(find(event_ts_round == time_idx_round))
            
            % add time_idx_round to event timeline:
            event_times_all=[event_times_all cur_time+time_idx/1000];
            event_numdetect_all=[event_numdetect_all num_detects];
            
            filenames_txt{length(filenames_txt)+1}= filename;
                % only show when number of detects > 2 (i.e. in at least 3 channels)
            if num_detects > min_select
                close all
                [b,a]=butter(6,450/(1000/2),'low');
                [c,d]=butter(6,5/(1000/2),'high');
                event_start=time_idx/FS-1;
                event_end=time_idx/FS+1;
                f=figure(4);
                f_ax=axes;
               % set(f, 'Visible', 'off')

                hold on
                g=figure(5);
                g_ax=axes;
              %  set(g, 'Visible', 'off')
                hold on
                samp_start = max(1,floor(event_start*FS));
                samp_end=min(length(data)-1,floor(event_end*FS));
                data_event = data(samp_start:samp_end,:);
                for mm=1:32

                    data_filt = filter(b,a,data(samp_start:samp_end,mm));%-data(:,32));
                    data_filt = filter(c,d,data_filt);
                    data_event(:,mm)=data_filt;
                   % data_filt=data_filt;
                    if mm<16
                        set(0, 'CurrentFigure', f)
                        subplot(10,3,mm)
                        specgram(data_filt./median(abs(data_filt)),100,FS,100,95)
                        title(['Ch = ' num2str(mm-1)])
                        caxis([0,40])
                        ylim([0,100])
                        %xlim([event_start,event_end])
                        set(0, 'CurrentFigure', g)
                        subplot(10,3,mm,'align')
                        plot((t(samp_start:samp_end)-t(1)),data_filt);
                        title(['Ch = ' num2str(mm-1)])
                        xlim([event_start,event_end])
                    elseif mm>16 && mm < 32
                        set(0, 'CurrentFigure', f)
                        subplot(10,3,mm-1)
                        specgram(data_filt./median(abs(data_filt)),100,FS,100,95)
                        title(['Ch = ' num2str(mm-1)])
                         caxis([0,40])
                        ylim([0,100])
                       % xlim([event_start,event_end])
                        set(0, 'CurrentFigure', g)
                        subplot(10,3,mm-1,'align')
                       plot((t(samp_start:samp_end)-t(1)),data_filt);
                        title(['Ch = ' num2str(mm-1)])
                        xlim([event_start,event_end])
                    end


                end
                
                % save event figures:
        set(0, 'CurrentFigure', f)
        set(f, 'Position', get(0, 'Screensize'));
            saveas(f,[mat_loc '/spec' directory(ii).name(1:end-4) 'plus' num2str(event_start) 's_' num2str(num_detects) '.png'])

            set(0, 'CurrentFigure', g)
            set(g, 'Position', get(0, 'Screensize'));
            saveas(g,[mat_loc '/ts' directory(ii).name(1:end-4) 'plus' num2str(event_start) 's_' num2str(num_detects) '.png'])
    % save event:
    t_event = t(samp_start:samp_end);
    save([mat_loc directory(ii).name(1:end-4) 'plus' num2str(event_start) 's_' num2str(num_detects) '.mat'],'t_event','data_event')
  % pause
    end
        end
% save over a .mat with all of the detections!

if floor(ii/100)==ceil(ii/100)
save([mat_loc '../eventcounts_' nameout '.mat'],'event_times_all','event_numdetect_all','filenames_txt')
end
        end
    end

end
%% 
save([mat_loc '../eventcounts_' nameout '.mat'],'event_times_all','event_numdetect_all','filenames_txt');

%% Trying with the .mats saved
