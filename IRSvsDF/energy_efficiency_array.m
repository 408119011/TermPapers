%% ECN 618 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006

function [EE_SISO,EE_DF,EE_IRS] = energy_efficiency_array(achievable_rate,B,efficiency_poweramp,power_dissipation_source, ...
    power_dissipation_destination,power_dissipation_relay,power_dissipation_perelement_irs,Noise_var,channelgain_sr,channelgain_rd,channelgain_sd,RefCoeff)
EE_SISO = zeros(length(achievable_rate),1);
EE_IRS  = zeros(length(achievable_rate),1);
EE_DF   = zeros(length(achievable_rate),1);
for i=1:length(achievable_rate)
    [EE_SISO(i),EE_DF(i),EE_IRS(i)] = energy_efficiency(achievable_rate(i),B,efficiency_poweramp,power_dissipation_source, ...
        power_dissipation_destination,power_dissipation_relay,power_dissipation_perelement_irs,Noise_var,channelgain_sr,channelgain_rd,channelgain_sd,RefCoeff);

end
end