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

% alternate path: 40s overlap?
% 10s overlap: '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_10s/N1/'
% On the BIGTOP

csv_dir1 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_40s/N1/';
csv_dir2 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_40s/N2/';
csv_dir3 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_40s/N3/';
csv_dir4 = '~/Downloads/SIDEX/data/box/Data_sigicom_2020/Processed_Data/results_40s/N4/';


%% On my laptop
%~/Downloads/SIDEx/data/Data_sigicom_2020/Processed_Data/results_40s/
csv_dir1 = '~/Downloads/SIDEx/data/Data_sigicom_2020/Processed_Data/results_40s/N1/';
csv_dir2 = '~/Downloads/SIDEx/data/Data_sigicom_2020/Processed_Data/results_40s/N2/';
csv_dir3 = '~/Downloads/SIDEx/data/Data_sigicom_2020/Processed_Data/results_40s/N3/';
csv_dir4 = '~/Downloads/SIDEx/data/Data_sigicom_2020/Processed_Data/results_40s/N4/';

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
     
     open_file=['csv_files' num2str(nodeindex)]      % change the name of the csv_file we're looking at
     event_name = strings(size(open_file)); % change csv_files# as neede                                         % just to check that its the right one
 end
 
 %% Change the number of the file for each run
 
 for i=1:size(csv_files1)
     titlefile = csv_files1(i).name;             % ith entry of event_name is a name from csv_files
     divide_title = split(titlefile,'_');       % divid the name by the '_' character
     event_name(i) = string(divide_title{1});       % only save the date of occurance, works!
     
 end
 
 % Doing the histc for each instance
 
 [edges,~,i] = unique(event_name);
 count = accumarray(i(:),1,[numel(edges),1]);
 
 % Change below to match csv files
 date1 = datetime(edges);
 instances1 = count;
 log_instances1= log10(instances1);
 
%% Process in plot_sigicom_timeseries
% save as a .mat, load here and plot over weather, will need to re run
% derivs if needed
% use load(10s_logN1-N4.mat) 

index1 = datefind(date1,weather_date);
index2 = datefind(date2,weather_date);
index3 = datefind(date3,weather_date);
index4 = datefind(date4,weather_date);

trigger1 =zeros(size(weather_date));
trigger2 =zeros(size(weather_date));
trigger3 =zeros(size(weather_date));
trigger4 =zeros(size(weather_date));

for ii = 1:size(index1)
   nn = index1(ii);
   trigger1(nn) = 1;
end

for ii = 1:size(index2)
   nn = index2(ii);
   trigger2(nn) = 2;
end

for ii = 1:size(index3)
   nn = index3(ii);
   trigger3(nn) = 3;
end

for ii = 1:size(index4)
   nn = index4(ii);
   trigger4(nn) = 4;
end



 
 %% Plot Each log Instance
 % one
 % yyaxis right
 
% subplot(2,1,2) % dont respecify since holding ON
 
xticks(weekticks)
 datetick('x', 'dd-mmm')
 %scatter(date1,instances1,70,'blue','o')
 scatter(date1,instances1,70,'blue','o')
 hold on

 % two
% scatter(date2,instances2,70,'blue','x')
scatter(date2,instances2,70,'red','x')
 hold on
 
 % three
% scatter(date3,instances3,70,'blue','s')
  scatter(date3,instances3,70,'green','s')
 hold on
 
 % four
% scatter(date4,instances4,70,'blue','^')
  scatter(date4,instances4,70,'black','^')
  ylabel('Number of Events Detected')

  set(gca, 'YScale', 'log')
  
  
 hold on

%% Three: Plot all the instances WITH temperature data
%% Weather data read in DONT RUN AGAIN

original_weather = readtable('SIDEx 2020 Weather.xlsx','Range','A2:D123');
% top line of each is the name of the

% Weather data processing
weather_date = original_weather.Jan; %(:,1); % first column is datetime
weather_maxC = (5/9)*(original_weather.Max - 32);  % second column is max
weather_minC = (5/9)*(original_weather.Min - 32); % fourth column is the MIN temp
weather_avgC = (5/9)*(original_weather.Avg - 32);  % third column is AVG use this one

%% creating datetime array for weeks
startdate = weather_date(1);
enddate =  weather_date(121);
weekticks = linspace(startdate,enddate,17);
%weekticks = datenum(weekticks);

%% WORKS: weather plotting, ONLY RUN THIS ONE

load('weather_data.mat');

%% Temperature derivatives

% empty array for holding
temp_derivs = double.empty(0,length(weather_avgC));
nn_add = 1;
% 

for nn=2:length(weather_avgC)   % for finding diff
    
     temp_derivs(nn_add) = weather_avgC(nn)-weather_avgC(nn-1);
     nn_add = nn_add+1;   
        
end
    
temp_derivs = transpose(temp_derivs);
weather_deriv = weather_date(1:120);

%% Works! Daily AVERAGE Temp

%yyaxis left
%plot(datenum(weather_date), weather_avgC);
%subplot(2,1,2)

plot(weather_date, weather_avgC);
hold off
%title('Weather and Number of Events Per Day')

xticks(weekticks)
datetick('x', 'dd-mmm','keepticks')

ylabel('Average Temperature (°C)')
xlabel('Dates, 2020')
xlim('auto');
%ax = gca;

yline(-10)
yline(-20)
legend('Avg. Temp °C','-10 °C', '-20 °C','1 Node','2 Node','3 Node','4 Node')  %,'Interpreter','latex')
title('Events Detected by Node and Temperature')



%% Plotting derivative
% make the right axis the dT/dt (difference in avg temp each day), also
% specify starting avg temp
% make the left axis the whole plot? with also nodes?
% x axis stays as the days it currently is

yyaxis left
plot(weather_deriv, temp_derivs);
datetick('x', 'dd-mmm-yyyy','keepticks')
ylabel('Change in Temperature (C)')
xlabel('Dates, 2020')
xlim('auto')
legend('1 Node','2 Node','3 Node','4 Node','$\frac{dT}{dt}$','Interpreter','latex')

title('Temperature Effect on Events (10s Overlap)')
set(gca, 'FontSize',15)
hold off

%% one big code of subplotting part ONE

subplot(2,1,1)
% need a count of # of triggers per day dived by number of nodes?h
scatter(weather_date,trigger1,70)
hold on
scatter(weather_date,trigger2,70)
hold on
scatter(weather_date,trigger3,70)
hold on
scatter(weather_date,trigger4,70)

xticks(weekticks)
datetick('x', 'dd-mmm','keepticks')
ylabel('Number of Nodes Triggered')
title('Nodes Triggered per Day, 10s Overlap')
xlabel('Dates, 2020')
ylim([0 5])
xtickangle(45)
set(gca, 'FontSize',12)
xlim('auto');


% one big code of subplotting part TWO

subplot(2,1,2) % dont respecify since holding ON

yyaxis right
datetick('x', 'dd-mmm')
%scatter(date1,instances1,70,'blue','o')
scatter(date1,instances1,70,'blue','o')
hold on

% two
% scatter(date2,instances2,70,'blue','x')
scatter(date2,instances2,70,'red','x')
hold on

% three
% scatter(date3,instances3,70,'blue','s')
scatter(date3,instances3,70,'green','s')
hold on

% four
% scatter(date4,instances4,70,'blue','^')
scatter(date4,instances4,70,'black','^')
ylabel('Number of Events Detected')
hold on

set(gca, 'YScale', 'log')

yyaxis left
plot(weather_date, weather_avgC);
hold off
%title('Weather and Number of Events Per Day')

xticks(weekticks)
datetick('x', 'dd-mmm','keepticks')
ylabel('Average Temperature (°C)')
xlabel('Dates, 2020')
xlim('auto');
xtickangle(45)
%ax = gca;

legend('Avg. Temp °C','1 Node','2 Node','3 Node','4 Node')  %,'Interpreter','latex')
title('Number Events Detected by Node and Weather')
set(gca, 'FontSize',12)
