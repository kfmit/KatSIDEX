close all
clear all
% Note: first import your csv as a matrix and call it "data"

% set directory:
% csv_dir = '~/data/sidex/sidex_2020_01_25/SigicomData/'
%% RUN if on bigtop
csv_dir = '~/Downloads/SIDEX/data/raw_transients/Data_sigicom_2020/Processed_Data/raw_transients/'

%% RUN if on laptop
csv_dir = '~/Downloads/SIDEx/data/Data_sigicom_2020/Processed_Data/raw_transients/'

%%

csv_files = dir([csv_dir '*.csv']);
standalone_GPS =[
     71.3357 -156.3982
     71.3299 -156.4016
   71.3345 -156.4165
   71.3333 -156.4081];
Node_IDS=[103212,103637,103636,103208];
Data_mat = zeros(1,4);
for ii=1:length(csv_files)
    % read in data:
    data_local=readtable([csv_files(ii).folder '/' csv_files(ii).name]);
   % data_local(cellfun('isempty',data_local)) = {nan}
    % parse variables names:
    Var_names = data_local.Properties.VariableNames;

    % prep the data matrix appropriately:
    mat_len = size(Data_mat,1)+1;
    local_len = length(data_local.offset)
    % original: local_len = length(data_local.ts)
    
    %Data_mat(mat_len:mat_len+local_len,1:(length(Node_IDS))*3+1)=nan;
    % put in time data:
    time_posix=posixtime(datetime(data_local.offset,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'));
    Data_mat(mat_len:mat_len+local_len-1,1)=time_posix(1:local_len);
%figure(2)
%hold on
    for kk=2:length(Var_names)
        % get Node ID, direction, time delay, timestamp
        % also add to node_IDS 
        vars=strsplit(Var_names{kk},'_');
        rowid=str2num(vars{2});
        direction = vars{end}
        
        tstamp = str2num(vars{1}(2:end)) + str2num(['0.' vars{end-2}]);
        if length(find(Node_IDS==rowid,1))<1
            % add to Node IDS
            Node_IDS=[Node_IDS rowid];
        end
        % insert data from column into appropriate part of Data_mat:
        node_idx = find(Node_IDS==rowid,1);
        
        if direction == 'V'
            mat_idx = (node_idx-1)*3 + 2;
        elseif direction == 'L'
            mat_idx = (node_idx-1)*3 + 3;
        elseif direction == 'T'
            mat_idx= (node_idx-1)*3+4;
        end
        mat_idx
        %%
        % insert into data matrix:
        if ~isa(table2array(data_local(1,kk)),'double')
            Data_mat(mat_len:mat_len+local_len-1,mat_idx) = str2double(table2array(data_local(1:local_len,kk)));
            
        else
            Data_mat(mat_len:mat_len+local_len-1,mat_idx) = table2array(data_local(1:local_len,kk));
        end
        %%
        %plot(table2array(data_local(1:local_len,kk)));
        %hold on
        %plot(Data_mat(mat_len:mat_len+local_len-1,mat_idx))
        %pause
    end
end
%%
% plot all channels:


%%
% load GPS data
%cal_file = '~/data/sidex/sidex_2020_01_25/GPSLOG00.TXT'
cal_file = '~/Downloads/SIDEX/data/GPS_Data/GPSLOG_2020-03-02.txt'


cal_fs=5;

A=readtable(cal_file);
cal_time=table2array(A(:,1));
cal_date = table2array((A(:,2)));
cal_lat=table2array(A(:,3));
cal_long=table2array(A(:,4));

a_x = table2array(A(:,6));
a_y=table2array(A(:,7));
a_z=table2array(A(:,8));
abs_accel = sqrt(a_x.^2+a_y.^2+a_z.^2);
cal_t =datetime('2020-01-25')+cal_time
cal_time_posix =posixtime(cal_t);
%%
% plot all 
figure(1)
subplot(2,1,1)
hold on
Data_mat(1,1)=Data_mat(2,1);
 for ii=2:(length(Node_IDS)*3+1)
     % plot the data:
     plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60,Data_mat(:,ii)/max(Data_mat(:,ii))+ii-2,'-')
     
     ylabel('Norm Channel + Voltage')
 end
 xlim([0,5000])
 subplot(2,1,2)
 hold on
 plot(cal_time_posix-cal_time_posix(1),abs_accel)
 xlabel('Time (s)')
 ylabel('Acceleration (m/s/s)')
 xlim([0,5000])

%%
% for each node: create a plot w/ L, T, V, acceleration, and range to jump
R_node(:,1) =getLocalRange(cal_lat,cal_long,standalone_GPS(1,1),standalone_GPS(1,2));
R_node(:,2) = getLocalRange(cal_lat,cal_long,standalone_GPS(2,1),standalone_GPS(2,2));
R_node(:,3)=getLocalRange(cal_lat,cal_long,standalone_GPS(3,1),standalone_GPS(3,2));
R_node(:,4) = getLocalRange(cal_lat,cal_long,standalone_GPS(4,1),standalone_GPS(4,2));
%%
min_lim=530;
max_lim=550;
for ll=4:length(Node_IDS)
    %figure('units','normalized','outerposition',[0 0 1 1])
    figure
    subplot(5,1,1)
    plot(cal_time_posix-cal_time_posix(1),abs_accel)
    xlim([min_lim,max_lim])
    title(['Node ID = ' num2str(Node_IDS(ll))])
    ylabel('Accel m/s^2')
    xlabel('Time (s)')
    subplot(5,1,2)
    plot(cal_time_posix-cal_time_posix(1),R_node(:,ll))
    xlim([min_lim,max_lim])
    title('Range (m)')
    ylabel('Range (m)')
    subplot(5,1,3)
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60,Data_mat(:,(ll-1)*3+2)/max(Data_mat(:,(ll-1)*3+2)),'-')
    
    ylabel('V Voltage')
    xlim([min_lim,max_lim])
    subplot(5,1,4)
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60,Data_mat(:,(ll-1)*3+3)/max(Data_mat(:,(ll-1)*3+3)),'-')
    xlim([min_lim,max_lim])
    ylabel('L Voltage')
    
    subplot(5,1,5)
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60,Data_mat(:,(ll-1)*3+4)/max(Data_mat(:,(ll-1)*3+4)),'-')
    xlim([min_lim,max_lim])
    ylabel('T Voltage')
    xlabel('Time (s)')
    
    %saveas(gcf,[csv_dir 'CalNode_' num2str(ll-1) '.png'])   
    %saveas(gcf,[csv_dir 'CalNode_' num2str(ll-1) '.fig'])   
end
%%
% read in and plot v. node 3 data from cabled array:
prefix = '~/data/sidex/sidex_2020_01_25/';
directory = dir([prefix 'Sidex*.txt']);
N=size(directory)


data = [];
FS = 1000;
stind=1;%length(directory)-10;
eind=length(directory);
startna=split(directory(stind).name,'Sidex_');
date_start=startna{2};

for i = stind:eind
    filename = [prefix directory(i).name]
    err=0;
    try 
        M = dlmread(filename, ',', 2, 0);
    catch
        disp('file error')
        M=[];
    end
    size(M)
   % pause
    if size(M,2)==16
    
    data = [data; M(1:(length(M)-1),:)];
    end
    %figure(2)
    %clf
    %hold on
    ch_plot=1:16;%[1,2,3,4,5,6,7,8,9,10,11,12,14,15]%1:16;%[1,2,3,4,5,6,10,11,12,13,15,16]
   % for ii=ch_plot
    %    
     %   plot(M(:,ii)+(ii-1)*max(abs(M(:,13))))
   % end
    %startna=split(directory(stind).name,'Sidex_');
    %date_start=startna{2};
end
% get start time/date, convert to epoch: (for now, look up on GPS file)
%%
t_start = 194038;
hour = floor(t_start/10000); 
minute = floor((t_start-hour*10000)/100); 
second = round(t_start-floor(t_start/100)*100); 
t_datetime = datetime(2020,1,25,hour,minute,second);
t_epoch=posixtime(t_datetime);
offset=-0.29;% offset of pps to second%419.75;
t_vec_cabled =  (0:1/FS:size(data,1)/FS - 1/FS) + t_epoch+offset;
figure
subplot(3,1,1)
hold on
    plot(t_vec_cabled-cal_time_posix(1)-300,data(:,10)/max(data(:,10)),'--');
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60-300,Data_mat(:,(ll-1)*3+2)/max(Data_mat(:,(ll-1)*3+2)),'-')
    
    xlim([20,30])
    ylabel('V Norm voltage')
    legend({'Cabled','Stand-alone'})
subplot(3,1,2)
hold on
plot(t_vec_cabled-cal_time_posix(1)-300,data(:,11)/max(data(:,12)),'--');
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60-300,Data_mat(:,(ll-1)*3+3)/max(Data_mat(:,(ll-1)*3+3)),'-')
    
    xlim([20,30])
    ylabel('L Norm voltage')
    subplot(3,1,3)
    hold on
    plot(t_vec_cabled-cal_time_posix(1)-300,data(:,12)/max(data(:,11)),'--');
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60-300,Data_mat(:,(ll-1)*3+4)/max(Data_mat(:,(ll-1)*3+4)),'-')
    
    xlim([20,30])
    ylabel('T Norm voltage')
    % plot pps:
    %plot(t_vec_cabled-cal_time_posix(1),data(:,13)/max(data(:,13))+3);
    %%
pause
% plot cabled v. not cabled:
figure;
ll=4;
    subplot(5,1,1)
    plot(cal_time_posix-cal_time_posix(1),abs_accel)
   % xlim([min_lim,max_lim])
    title(['Node ID = ' num2str(Node_IDS(4))])
    ylabel('Accel m/s^2')
    xlabel('Time (s)')
    subplot(5,1,2)
    plot(cal_time_posix-cal_time_posix(1),R_node(:,ll))
    %xlim([min_lim,max_lim])
    title('Range (m)')
    ylabel('Range (m)')
    
    subplot(5,1,3)
    hold on
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60,Data_mat(:,(ll-1)*3+2)/max(Data_mat(:,(ll-1)*3+2)),'-')
    plot(t_vec_cabled-cal_time_posix(1),data(:,10)/max(data(:,10)));
    ylabel('V Voltage')
    %xlim([min_lim,max_lim])
    subplot(5,1,4)
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60,Data_mat(:,(ll-1)*3+3)/max(Data_mat(:,(ll-1)*3+3)),'-')
    hold on
    plot(t_vec_cabled-cal_time_posix(1),data(:,11)/max(data(:,11)))
    xlim([min_lim,max_lim])
    ylabel('L Voltage')
    
    subplot(5,1,5)
    plot(Data_mat(:,1)-cal_time_posix(1)+5*60*60,Data_mat(:,(ll-1)*3+4)/max(Data_mat(:,(ll-1)*3+4)),'-')
    xlim([min_lim,max_lim])
    ylabel('T Voltage')
    xlabel('Time (s)')
    hold on
    plot(t_vec_cabled-cal_time_posix(1),data(:,12)/max(data(:,12)))

%%
save([csv_dir 'caldata_sigicom.mat'])