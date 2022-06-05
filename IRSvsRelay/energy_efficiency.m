%% ECN 618 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006

function [EE_SISO,EE_DF,EE_IRS] = energy_efficiency(achievable_rate,B,efficiency_poweramp,power_dissipation_source,power_dissipation_destination, ... 
                power_dissipation_relay,power_dissipation_perelement_irs,Noise_var,channelgain_sr,channelgain_rd,channelgain_sd,RefCoeff)
% energy efficiency calculation through capacity
sinr          = 2^(achievable_rate)-1;           % SISO and IRS
sinr_df       = 2^(2*achievable_rate)-1;         % DF relaying
power_siso    = sinr*Noise_var/channelgain_sd;   % SISO power

% realy power
power_df = sinr_df*Noise_var*(channelgain_sr+channelgain_rd-channelgain_sd)/(2*channelgain_rd*channelgain_sr);

% optimal number of elements needed in irs power calculation
Nopt = (2*sinr*Noise_var/(RefCoeff^2*channelgain_sr*channelgain_rd*power_dissipation_perelement_irs)).^(1/3) ...
                - sqrt(channelgain_sd/(channelgain_sr*channelgain_rd))/RefCoeff;

% optoimised energy efficiency in siso 
EE_SISO = 1000*B*achievable_rate/(power_siso/efficiency_poweramp + power_dissipation_source + power_dissipation_destination);

% optimisednenergy efficincy in relay
EE_DF = 1000*B*achievable_rate/(power_df/efficiency_poweramp + power_dissipation_source/2 ...
                                      + power_dissipation_destination + power_dissipation_relay);
if Nopt<0
    Nopt = 0;
end

% optimised power in irs
power_irs = sinr*Noise_var/(sqrt(channelgain_sd) + Nopt*RefCoeff*sqrt(channelgain_sr*channelgain_rd))^2;

% ee in irs
EE_IRS = 1000*B*achievable_rate/(power_irs/efficiency_poweramp + power_dissipation_source ... 
                     + power_dissipation_destination + Nopt*power_dissipation_perelement_irs);

end