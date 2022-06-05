%% ECN 620 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006 

function [SNR, theta, h, ratio1, ratio2, SNRPlot] = Proposed_Algorithm_multilayer(g, f, Sigma, Pmax, rpt, R, attenuation, N1, EUR_threshold)
N = size(g, 1);          % number of elements
K = size(f, 2)-(R-1)*N;  % users
v = ones(size(g, 2), 1); % recombining 
theta = exp(1j*2*pi*rand(N, R));% randomely generated TPS matrix
SNRPlot = zeros(rpt, 1);        % snr vector
%Algorithm implementation
% snr calculation as given in paper through the proposed algorithm
for t = 1:rpt                   % number of point to be snr calculted on 
    a = g*v;
    for r = R:-1:2
        a = f(:, K+1+N*(r-2):K+N*(r-1))'*diag(theta(:, r))'*a; %optimization exprerssion of paper
    end
    a = f(:, 1:K)'*diag(theta(:, 1))'*a;
    w = a/norm(a);              % unit power beam forming vector
    b = v'*g';
    for r = R:-1:2
        b = b*diag(theta(:, r))*f(:, K+1+N*(r-2):K+N*(r-1));
    end
    b = b*diag(f(:, 1:K)*w);
    theta(:, 1) = exp(1j*angle(b')); % tps optimization
    %tps optimization
    for l = 2:R
        b = v'*g';
        for r = R:-1:l+1
            b = b*diag(theta(:, r))*f(:, K+1+N*(r-2):K+N*(r-1));
        end
        DiagMat = f(:, K+1+N*(l-2):K+N*(l-1));
        for r = l-1:-1:2
            DiagMat = DiagMat*diag(theta(:, r))*f(:, K+1+N*(r-2):K+N*(r-1));
        end
        DiagMat = DiagMat*diag(theta(:, 1))*f(:, 1:K)*w;
        b = b*diag(DiagMat);
        theta(:, l) = exp(1j*angle(b')); % optimised tps
    end
    Prd = diag(theta(:, 1))*f(:, 1:K); % recombinig optimization
    % recombining optimization
    for r = R:-1:2
        Prd = diag(theta(:, r))*f(:, K+1+N*(r-2):K+N*(r-1))*Prd;
    end
    [v, ~] = eigs(g'*Prd*w*w'*Prd'*g, 1); %optimised recombining vector
    SNR = abs(v'*g'*Prd*w)^2/norm(v)^2/Sigma;
    SNRPlot(t) = SNR;
end
SNR = SNR*Pmax;

Width1 = N1(1);
Width2 = N1(2);
h = zeros(Width1, Width2, R);
PlotPrd = f(:, 1:K)*w;
h(:, :, 1) = reshape(PlotPrd, Width2, Width1)';
for r = 2:R
    PlotPrd = f(:, K+1+N*(r-2):K+N*(r-1))*diag(theta(:, r-1))*PlotPrd;
    PlotPrd = PlotPrd*attenuation;
    h(:, :, r) = reshape(PlotPrd, Width2, Width1)';
end

%% energy plot on both of the layers
%%% plot in R figures

for r = 1:R
    figure;
    imagesc(pow2db(abs(h(:, :, r)).^2), [-40 -10]);
    colormap(jet);
    axis equal;
    axis off;
    title(['Fig. 7. UC-RIS: Layer ', num2str(r)], 'Fontsize', 15);
end

Energy1 = sum(abs(h(:, :, 1)).^2, 'all');
Energy2 = sum(abs(h(:, :, 2)).^2, 'all');

ratio1 = sum(abs(h(:, :, 1)).^2>Energy1/numel(h(:, :, 1))*EUR_threshold, 'all')/numel(h(:, :, 1));
ratio2 = sum(abs(h(:, :, 2)).^2>Energy2/numel(h(:, :, 2))*EUR_threshold, 'all')/numel(h(:, :, 2));
end