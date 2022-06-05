%% ECN 618 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006

function output = pathloss_NLOS(distance,fc)
output = 10.^(0.1*(-22.7-26*log10(fc)-36.7*log10(distance)));
end