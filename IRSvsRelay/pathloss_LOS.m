%% ECN 618 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006

function output = pathloss_LOS(distance,fc)
output = 10.^(0.1*(-28-20*log10(fc)-22*log10(distance)));
end