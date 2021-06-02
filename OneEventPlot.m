%% Using just the V for now, copy for T and L or jsut change the 'v' in each file
clear
clc
close all

% this could also be accomplished with a loop

tranv3 = readtable('transient_04082020_mp_3_c22_103208_v_2020-01-26-070006146.xlsx', 'Range', 'A5:B168964'); %read in the v from node 3
tranv4 = readtable('transient_04082020_mp_4_c22_103212_v_2020-01-26-070005106.xlsx', 'Range', 'A5:B168964'); %read in the v from node 4
tranv5 = readtable('transient_04082020_mp_5_c22_103636_v_2020-01-26-070005156.xlsx', 'Range', 'A5:B164868'); % read in the v from node 5
tranv6 = readtable('transient_04082020_mp_6_c22_103637_v_2020-01-26-070005583.xlsx', 'Range', 'A5:B168964'); % read in v from node 6

%% Assigning Variables. All of this could be done with a for loop but there were only four so i went with less troubleshooting

timestamp3 = datetime(2020,01,26,07,00,61,46); % literally written in, matlab only reads cell as NaN
timev3 = timestamp3 + seconds(table2array(tranv3(:,1))); %create time x using column 1 for L
velv3 = table2array(tranv3(:,2)); %create l using column 2 for L

timestamp4 = datetime(2020,01,26,07,00,51,06);
timev4 = timestamp4 + seconds(table2array(tranv4(:,1))); % same with node 4
velv4 = table2array(tranv4(:,2));

timestamp6 = datetime(2020,01,26,07,00,51,56);
timev6 = timestamp6 + seconds(table2array(tranv6(:,1))); % same with node 6
velv6 = table2array(tranv6(:,2));

timestamp5 = datetime(2020,01,26,07,00,55,83);
timev5 = timestamp5 + seconds(table2array(tranv5(:,1))); % same with node 5
velv5 = table2array(tranv5(:,2));


%% Plotting the graphs
figure
tiledlayout(4,1)

nexttile
plot(timev3,velv3)
title('Node 3')
xlabel('Time, (s)')
ylabel('Velocity (in/s)')
xlim([timestamp3 timev5(164864)])

nexttile
plot(timev4,velv4)
title('Node 4')
xlabel('Time, (s)')
ylabel('Velocity (in/s)')
xlim([timestamp3 timev5(164864)])

nexttile
plot(timev6,velv6)
title('Node 6')
xlabel('Time, (s)')
ylabel('Velocity (in/s)')
xlim([timestamp3 timev5(164864)])

nexttile
plot(timev5,velv5)
title('Node 5')
xlabel('Time, (s)')
ylabel('Velocity (in/s)')
xlim([timestamp3 timev5(164864)])



%% Secondary Graph with One continuous plot
figure(2)
hold on
plot(timev3,velv3)
plot(timev4,velv4)
plot(timev5,velv5)
plot(timev6,velv6)
title('Combined Nodes')
xlabel('Time, (s)')
ylabel('Velocity (in/s)')
xlim([timestamp3 timev5(164864)])
legend('Node 3', 'Node 4', 'Node 6', 'Node 5')
