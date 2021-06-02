%% 02JUN2021
% Kat Fung
% Time Series

clc
clear
close all

% One: Import all of the csv from Sigicom2020
% At least their names, bc you can pull the timestamp from there
csv_dir = '~/Downloads/SIDEX/data/raw_transients/Data_sigicom_2020/2510a436-69d3-4757-8f35-119bf88f566f/raw_transients/'

csv_files = dir([csv_dir '*.csv']);
standalone_GPS =[
     71.3357 -156.3982
     71.3299 -156.4016
   71.3345 -156.4165
   71.3333 -156.4081];
Node_IDS=[103212,103637,103636,103208];
Data_mat = zeros(1,4);

% Note to self: may need to make a for loop for adding the offset onto the
% epoch time to fix the errors

%%
for ii=1:length(csv_files)
    % read in data:
    data_local=readtable([csv_files(ii).folder '/' csv_files(ii).name]);
   % data_local(cellfun('isempty',data_local)) = {nan}
    % parse variables names:
    Var_names = data_local.Properties.VariableNames;

    % prep the data matrix appropriately:
    mat_len = size(Data_mat,1)+1;
    % WHat the living hell is OFFSET? possible the difference from the
    % timestamp?
    local_len = length(data_local.offset)
    % orig: local_len = length(data_local.ts)
    
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

%% Two: Plot all of these instances of triggering on a timeseries plot

% things I will probably want: x axis is date in day-mon-yr form, left y is
% # of events on that day, right y is temperature maybe?

%% Three: Plot all the instances WITH temperature data

% pull the temp data that you used with the old graph
