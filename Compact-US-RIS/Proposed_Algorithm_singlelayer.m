%% ECN 620 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006 

function [SNR3, SNR5, h, theta, SNR3_Plot, SNR5_Plot, ratio]  = Proposed_Algorithm_singlelayer(g, f, Sigma, Pmax, rpt, attenuation, N1, EUR_threshold)
v = ones(size(g, 2), 1);
theta = ones(size(g, 1), 1);
SNR3_Plot = zeros(rpt, 1);
SNR5_Plot = zeros(rpt, 1);
for t = 1:rpt
    a = f'*diag(theta)'*g*v;
    w = a/norm(a);
    b = v'*g'*diag(f*w);
    theta = exp(1j*angle(b'));
    [v, ~] = eigs(g'*diag(theta)*f*w*w'*f'*diag(theta)'*g, 1);
    SNR3 = abs(v'*g'*diag(theta)*f*w)^2/norm(v)^2/Sigma;
    SNR5 = abs(attenuation*v'*g'*diag(theta)*f*w)^2/norm(v)^2/Sigma;
    SNR3_Plot(t) = SNR3;
    SNR5_Plot(t) = SNR5;
end
SNR3 = SNR3*Pmax;
SNR5 = SNR5*Pmax;

Width1 = N1(1);
Width2 = N1(2);
h = reshape(f*w, Width2, Width1)';
figure;
imagesc(pow2db(abs(h).^2), [-40 -10]);
axis equal;
axis off;
title('Fig. 7.')
colormap(jet);
colorbar('Ticks', [-40 -30 -20 -10 0], 'FontSize', 20,'TickLabelInterpreter', 'latex');

Energy = sum(abs(h).^2, 'all');
ratio = sum(abs(h).^2>Energy/numel(h)*EUR_threshold, 'all')/numel(h);
end