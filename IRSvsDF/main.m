close all;
clear;
clc;

%% ECN 618 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006 

%% Hardware & channel Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 3GHz
% Bandwidth = 10 MHz
% figure of noise = 10 dB
% Reflection coefficient = 1
% Antenna gains = 5, 5, 0 dB in source, IRS, Receiver respectively
% efficiency of power amplifier = 0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc             =  3;                % carrier frequency in GHz
B              =  10e6;             % bandwidth in Hz
Noise_dB       =  10;               % noise figure in dB
RefCoeff       =  1;                % reflection coefficient

Gs = 10^(0.1*5);    % antenna gain at the Source,
Gr = 10^(0.1*5);    % antenna gain at the IRS,
Gd = 10^(0.1*0);    % antenna gain at the Destination

Noise_var_dBm         = -174 + 10*log10(B) + Noise_dB;     % noise power in dBm
Noise_var             =  10.^(0.1*Noise_var_dBm);          % linear noise power .
efficiency_poweramp   =  0.5;

%% Locations of components (In meters)
dist_source_irsnrelay             = 80; % distance between source and irs or relay as irs and relay is at same location
normal_dist_lineofSOURCEnIRS_DEST = 10; % perpendicular deisatnce between source and irs
parallel_dist_SOURCE_DEST         = 70; % parallel distance between source and destination

% distance between source and destination
dist_source_dest    =  distance_calc(normal_dist_lineofSOURCEnIRS_DEST,parallel_dist_SOURCE_DEST);% distance calculation between two points

% distance between relay and destination
dist_irsnrelay_dest =  distance_calc(normal_dist_lineofSOURCEnIRS_DEST,dist_source_irsnrelay-parallel_dist_SOURCE_DEST);

%% channels gains vs distance PLOT 1
% The channel gains are modeled using the 3GPP Urban Micro (UMi)
% from [17, Table B.1.2.1-1] with a carrier frequency of 3 GHz.

% distanec vector
d = 10:0.01:100;
channelgain_LOS  = 10*log10(pathloss_LOS(fc,d)*Gs*Gr);% channel gain los
channelgain_NLOS = 10*log10(pathloss_NLOS(fc,d)*Gs*Gr);% channel gain nlos

%% Plot
%channel gain plots
figure('Name','Typical channel gains as a function of the distance');
plot(d,channelgain_LOS,'k--','LineWidth',2);
hold on;
plot(d,channelgain_NLOS,'r--','LineWidth',2);
grid on;
title('Typical channel gains as a function of the distance');
xlabel('Distance ');
ylabel('Channel gain in dB');
legend('UMi-LOS','UMi-NLOS');


%% PLOT 2

%% channel gains
channelgain_sr      =  pathloss_LOS(dist_source_irsnrelay,fc)*Gs*Gr; % channel between source and relay
channelgain_rd      =  pathloss_LOS(dist_irsnrelay_dest,fc)*Gr*Gd;   % channel gain between source and destinatino
channelgain_sd      =  pathloss_NLOS(dist_source_dest,fc)*Gs*Gd;     % channel gain betwenn source and destination

power_dissipation_source          = 100;                             % (mW)% power dissipation in souce
power_dissipation_destination     = 100;                             % (mW)% in destination
power_dissipation_per_element_irs = 5;                               % (mW)% per element of irs
power_dissipation_relay           = 100;                             % (mW)% in relay


Nrange                    = [25 50 100 150];                         % range of elements in irs
achievable_rate           = [4 6];                                   % capacity constraint as given in paper
parallel_dist_SOURCE_DEST = 40:100;                                  % varied distance between source and destination

% Nmin = zeros(length(achievable_rater),1);

%% Transmit power
for i = 1:length(achievable_rate)
    %transmit power calculation for each capacity
    [TxP_IRS,TxP_DF,TxP_SISO,Nmin] = transmit_power(achievable_rate(i),parallel_dist_SOURCE_DEST, ...
                       dist_source_irsnrelay,normal_dist_lineofSOURCEnIRS_DEST,fc,Gs,Gr,Gd,Noise_var,Nrange,RefCoeff );
    
    figure('Name',' Transmit Power vs Disatance in all cases  ');
    hold on;
    plot(parallel_dist_SOURCE_DEST,10*log10(TxP_SISO),'k--','LineWidth',2);
    plot(parallel_dist_SOURCE_DEST,10*log10(TxP_DF),'b-.','LineWidth',2);
    for n = 1:length(Nrange)
        plot(parallel_dist_SOURCE_DEST,10*log10(TxP_IRS(:,n)),'r-','LineWidth',2);
    end
    str = sprintf(' %d bits/sec/Hz rate constraint', achievable_rate(i));
    title(str);
    xlabel('Distance d_1 [m] (parallel distance between source and destination)');
    ylabel('Transmit power in dBm');
    legend('SISO','IRS','DF relay');
    fprintf('Min no. of IRS ele. req to beat relaying with rate %f should be >= %d \n',achievable_rate(i),Nmin);
    
end

%% PLOT 3
%% Energy Efficiency

achievable_rate = [0.01 0.1:0.1:10]; % capacity vector

% energy efficiency array calculated
[EE_SISO,EE_DF,EE_IRS] = energy_efficiency_array(achievable_rate,B,efficiency_poweramp,power_dissipation_source,power_dissipation_destination, ...
               power_dissipation_relay,power_dissipation_per_element_irs,Noise_var,channelgain_sr,channelgain_rd,channelgain_sd,RefCoeff);

figure('Name','Energy Efficiency');
hold on;
plot(achievable_rate,EE_DF/1e6,'b-.','LineWidth',2);
plot(achievable_rate,EE_IRS/1e6,'r-','LineWidth',2);
plot(achievable_rate,EE_SISO/1e6,'k--','LineWidth',2);
title('Energy Efficiency');
xlabel('Achievable rate in bit/s/Hz');
ylabel('Energy efficiency Mbit/Joule');
legend('DF relay','IRS','SISO');
