%% ECN 618 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006

function  [P_IRS,P_DF,P_SISO,Nmin] = transmit_power(achievable_rate,distance,dist_source_irsnrelay,normal_dist_lineofSOURCEnIRS_DEST,fc,Gs,Gr,Gd,Noise_var,Nrange,RefCoeff )

% transmit power calculation by using capacity formula in all cases
P_IRS  = zeros(length(distance),length(Nrange));         % transmit power vector for irs
P_DF   = zeros(length(distance),1);                      % transmit power vector for relay
P_SISO = zeros(length(distance),1);                      % transmit power vector for siso

SINR    = 2^(achievable_rate)-1;                         % SISO and IRS
SINR_DF = 2^(2*achievable_rate)-1;                       % DF relaying

for i = 1:length(distance)

    %Compute distance between the source and destination
    d_SD = distance_calc(distance(i),normal_dist_lineofSOURCEnIRS_DEST);

    %Compute distance between the IRS/relay and destination
    d_RD = distance_calc(distance(i)-dist_source_irsnrelay,normal_dist_lineofSOURCEnIRS_DEST);

    channelgain_sr = pathloss_LOS(dist_source_irsnrelay,fc)*Gs*Gr;    % channel model betweeen source and relay
    channelgain_rd = pathloss_LOS(d_RD,fc)*Gr*Gd;                     % channel model betweeen relay and destinatino
    channelgain_sd = pathloss_NLOS(d_SD,fc)*Gs*Gd;                    % channel model betweeen source and destinatiion

    P_SISO(i) = SINR*Noise_var/channelgain_sd;                        % power of siso

    % power in irs
    P_IRS(i,:) = SINR*Noise_var./(sqrt(channelgain_sd) + Nrange*RefCoeff *sqrt(channelgain_sr*channelgain_rd)).^2;

    % power in relay
    if channelgain_sr>=channelgain_sd
        P_DF(i) = SINR_DF*Noise_var*(channelgain_sr+channelgain_rd-channelgain_sd)/(2*channelgain_rd*channelgain_sr);
    else
        P_DF(i) = SINR_DF*Noise_var/channelgain_sd;
    end
    if channelgain_sr>=channelgain_sd

        rho  = P_DF(i)/Noise_var;

        Nmin = sqrt((sqrt(1+rho*2*channelgain_rd*channelgain_sr/(channelgain_sr+channelgain_rd-channelgain_sd))-1)/(rho*channelgain_sr*channelgain_rd))-sqrt(channelgain_sd/(channelgain_sr*channelgain_rd));
    else
        Nmin = 1;
    end

end
end
