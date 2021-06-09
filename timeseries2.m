%% 02JUN2021
% Kat Fung
% Time Series

clc
clear
close all

%% One: Import all of the csv from Sigicom2020
% At least their names, bc you can pull the timestamp from there
% There are four folders within results_10 (N1-N4)
%csv_dir = '~/Downloads/SIDEX/data/raw_transients/Data_sigicom_2020/Processed_Data/results_10s/N1/'

csv_dir1 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_10s/N1/';
csv_dir2 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_10s/N2/';
csv_dir3 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_10s/N3/';
csv_dir4 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_10s/N4/';

%% DO NOT RUN AGAIN, just load the .mat of the files
%for nodenum=1:4
    
    %csv_dir=['csv_dir' num2str(nodenum)];
    %csv_name=['csv_files' num2str(nodenum)];
    csv_files1 = dir([csv_dir1 '*.csv']);
    csv_files2 = dir([csv_dir2 '*.csv']);
    csv_files3 = dir([csv_dir3 '*.csv']);
    csv_files4 = dir([csv_dir4 '*.csv']);
    
%end


standalone_GPS =[
     71.3357 -156.3982
     71.3299 -156.4016
     71.3345 -156.4165
     71.3333 -156.4081];
Node_IDS=[103212,103637,103636,103208];
Data_mat = zeros(1,4);

% Note to self: may need to make an edit inside for loop for adding the offset onto the
% epoch time to fix the errors

%% Use loop to read in all data from files, Not needed for now
% Expand if needed
for ii=1:length(csv_files)
    % read in data:
    data_local=readtable([csv_files(ii).folder '/' csv_files(ii).name]);
   % data_local(cellfun('isempty',data_local)) = {nan}
    % parse variables names:
    Var_names = data_local.Properties.VariableNames;

    % prep the data matrix appropriately:
    mat_len = size(Data_mat,1)+1;
    
    % STEPS TO FIX
    % 1) cut off the first part of every csv time
    % 2) add the offset to the epoch time
    
    local_len = length(data_local.ts)
    %local_len = length(data_local.offset)
    % orig: local_len = length(data_local.ts)
    
    %Data_mat(mat_len:mat_len+local_len,1:(length(Node_IDS))*3+1)=nan;
    % put in time data:
    time_posix=posixtime(datetime(data_local.ts,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'));
    Data_mat(mat_len:mat_len+local_len-1,1)=time_posix(1:local_len);
%figure(2)
%hold on

%% Loop to assign name values (doesn't work yet)
    for kk=2:length(Var_names)                          % for the length of var_names
        % get Node ID, direction, time delay, timestamp
        % also add to node_IDS 
        vars=strsplit(Var_names{kk},'_');
        rowid=str2num(vars{2});
        direction = vars{end}
       
        tstamp = str2num(vars{1}(2:end)) + str2num(['0.' vars{end-2}]);
        if length(find(Node_IDS==rowid,1))<1
            add to Node IDS
            Node_IDS=[Node_IDS rowid];
        end
        insert data from column into appropriate part of Data_mat:
        node_idx = find(Node_IDS==rowid,1);
        
        if direction == 'V'
            mat_idx = (node_idx-1)*3 + 2;
        elseif direction == 'L'
            mat_idx = (node_idx-1)*3 + 3;
        elseif direction == 'T'
            mat_idx= (node_idx-1)*3+4;
        end
        mat_idx
        %
        insert into data matrix:
        if ~isa(table2array(data_local(1,kk)),'double')
            Data_mat(mat_len:mat_len+local_len-1,mat_idx) = str2double(table2array(data_local(1:local_len,kk)));
            
        else
            Data_mat(mat_len:mat_len+local_len-1,mat_idx) = table2array(data_local(1:local_len,kk));
        end
        
    end
end

%% Part 1.1: use the "csv_dir" names to establish points (datalocal has all the vlt for later)
% Pull names

%event_name=char.empty(size(csv_files1),1);        % create empty array

for nodeindex=[1:4]
    
    csv_files=['csv_files' num2str(nodeindex)]      % change the name of the csv_file we're looking at
    pause;                                         % just to check that its the right one
    
   %% internal for loop to pick apart names
    for i=1:size(csv_files)
        titlefile = csv_files(i).name;             % ith entry of event_name is a name from csv_files
        divide_title = split(titlefile,'_');       % divid the name by the '_' character
        event_name(i) = divide_title{1};       % only save the date of occurance
        
    end
    
end
%% Hardcode way of doing above
event_name = strings(size(csv_files1)); % change csv_files# as needed

    for i=1:size(csv_files1)
        titlefile = csv_files1(i).name;             % ith entry of event_name is a name from csv_files
        divide_title = split(titlefile,'_');       % divid the name by the '_' character
        event_name(i) = string(divide_title{1});       % only save the date of occurance, works!
        
    end
    
%% Doing the histc for each instance

[edges,~,i] = unique(event_name);
date1 = datetime(edges);
count = accumarray(i(:),1,[numel(edges),1]);
instances1 = count;


%% Two: Plot all of these instances of triggering on a timeseries plot

% things I will probably want: x axis is date in day-mon-yr form, left y is
% # of events on that day, right y is temperature maybe?

%% ONE: plot the N1
yyaxis left
plot(date1,instances1)
datetick('x', 'dd-mmm-yyyy')
title('Weather & Number of Events Per Day')
hold on

%% TWO: plot the N2
yyaxis left
plot(date2,instances2)
datetick('x', 'dd-mmm-yyyy')
hold on

%% THREE: plot the N3
yyaxis left
plot(date3,instances3)
datetick('x', 'dd-mmm-yyyy')
hold on

%% FOUR: plot the N4
yyaxis left
plot(date3,instances3)
datetick('x', 'dd-mmm-yyyy')
hold on

%% Not a hardcode way of doing it
for nodeindex=1:4
    
    date=['date' num2str(nodeindex)];
    instance=['instance' num2str(nodeindex)];
    plot(date,instance)
    hold on
    
end


%% Three: Plot all the instances WITH temperature data
%% Weather data read in DONT RUN AGAIN

original_weather = readtable('SIDEx 2020 Weather.xlsx','Range','A2:D123');
% top line of each is the name of the

% Weather data processing
weather_date = original_weather.Jan; %(:,1); % first column is datetime
weather_maxC = (5/9)*(original_weather.Max - 32);  % second column is max
weather_minC = (5/9)*(original_weather.Min - 32); % fourth column is the MIN temp
weather_avgC = (5/9)*(original_weather.Avg - 32);  % third column is AVG use this one

%% WORKS: weather plotting, ONLY RUN THIS ONE

load('weather_data.mat');

%% Works!
yyaxis right 
%plot(datenum(weather_date), weather_avgC);
plot(weather_date, weather_avgC);
%title('Weather and Number of Events Per Day')
datetick('x', 'dd-mmm-yyyy')
ylabel('Average Temperature (C)')
xlabel('Dates, 2020')
xlim('auto')
legend('1 Node','2 Node','3 Node','4 Node','Avg Temp')
hold off