%% Kat Fung
%% Sigicom Sidex Data Processing Part 2
%%

clear all
close all
clc

% load the already collated data
load('original_caldata_sigicom.mat'); % contains variables used in original collation
load('Data_mat1.mat'); %contains ALL of the actual data in a mat from plot_sigicom_git statcaldata (big file)
load('dates collected.mat'); % has the dates I pulled
load('datesort.mat'); % has the uniques and num of events
%% Weather data read in DONT RUN AGAIN

original_weather = readtable('SIDEx 2020 Weather.xlsx','Range','A1:D123');

%% Weather data processing
weather_date = original_weather.Time; %(:,1); % first column is datetime
weather_max = (5/9)*(original_weather.Temperature__F_ - 32); % second column is max
weather_avg = (5/9)*(original_weather.Var3 -32); % third column is AVG use this one
weather_min = original_weather.Var4; % fourth column is the MIN temp

%% weather plotting
%% CONVERT TO CELSIUS

yyaxis right 
plot(weather_date,weather_avg);
title('Weather and Number of Events Per Day')
ylabel('Average Temperature (C)')
xlabel('Dates, 2020') 
hold on
%% Doing # of events
% column 1: utc epoch time of event
% column 2-4: VLT of Node 1
% column 5-7: VLT of Node 2
% column 8-10: VLT of Node 3
% column 11-13: VLT of Node 4

%% Plot the events and numbers
% Variables are from TestPlot1.m: num_events, ymd_mat, ymd_unique
unique_dates = datetime(ymd_unique);
 %make sure datetime converted sucessfully
% actual plot 
yyaxis left
scatter(unique_dates,num_events)
ylabel('Number of Events per Day')
xlim([min(unique_dates)-3 max(unique_dates)+3])

