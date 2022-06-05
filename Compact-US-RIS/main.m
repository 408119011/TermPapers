%% ECN 620 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006 

close all; 
clear;
clc;

%% Parameters
f = 2.5e9;                 % carrier frequency
c = 3e8;                   % velocity of light
lambda = c/f;              % wavelenght of wave
rpt = 12;                  % number of points to ploted in snr curve
Dis_BS2RIS = 20;           % distance between base station plane and ris layer plane
Dis_RIS2User = 0.02;       % perpendicular distance between RIS and user plane
Dis_Layer = 0.02;          % distance between layers 
M = 8;                     % number of base station antennas
K = 2;                     % number of user equipment antennas
R = 2;                     % number of layers
N1 = [12, 16];             % single IRS layer case elements
N2 = [8, 12];              % two irs layer case elements
Sigma = 1e-6;              % power of noise
loss = 0.8;                % energy loss of transmissive RIS
EUR_threshold = 1/6;       % power treshold in calculating EUR
RISWidth = 2;              % sidewidth of RIS element = lambda/RISWidth
  
Pmax_dB = -3:1:3;          % power in db
Pmax = db2pow(Pmax_dB);    % maximum power


%% Multi-Layer Case or two layer  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generation of positions of Antennas of user, base station and element Positions of two layer IRS  
  
%% 1.USER
  User_y = 0;             % y coordinate of user antenna
  User_Pos = zeros(K, 2); % (x,y) cordinate of K users
  % generating (x,y) for all user antennas
  for k = 1:K % user antenna index
    User_Pos(k, :) = [-(ceil(K/2)-1)*lambda/2-lambda/4+(k-1)*lambda/2, 0];% expression depend on index and lambda as given in paper
  end

%% 2.RIS position generation
RIS_y = Dis_RIS2User+Dis_Layer*(0:R-1);  % y coordinate of all RIS layers
RIS_Pos = zeros(N2(1)*N2(2), 2);         % (x,y) coordinate of all RIS layes
sideWidth = lambda/RISWidth;             % diastance between layers
% generating position of each elements in multi  layer
for n = 1:N2(1)*N2(2) % element index
    % X coordinate = -(N2(2)/2-1)*sideWidth-sideWidth/2+(n-(y-1)*N2(2)-1)*sideWidth
    % Y coordinate = (N2(1)/2-1)*sideWidth+sideWidth/2-(ceil(n/N2(2))-1)*sideWidth,
    % Y coordinate is independent of element index
    y = ceil(n/N2(2));
    x = n-(y-1)*N2(2);
    RIS_Pos(n, :) = [-(N2(2)/2-1)*sideWidth-sideWidth/2+(x-1)*sideWidth,(N2(1)/2-1)*sideWidth+sideWidth/2-(y-1)*sideWidth];% depend on index and lambda given in paper
end

%% 3.Base station position
BS_y = Dis_RIS2User*R+Dis_BS2RIS;   % y coordinate of base station antennas (is same)
BS_Pos = zeros(M, 2);               % (x,y) coordinate of base station antennas
for k = 1:M                         % base station antenna index
    BS_Pos(k, :) = [-(ceil(M/2)-1)*lambda/2-lambda/4+(k-1)*lambda/2, 0];% y coordinate is fixed, and x depend on lambda
end


%% Generate Channels

N = N2(1)*N2(2);% total number of elements in each layer

g = zeros(N, M);% channel between 1st layer of IRS to base station
%% 1.Channel matrix g between IRS and Base station
for n = 1:N     % element index of irs
    for m = 1:M % base station antenna index
        Dis = sqrt((RIS_Pos(n, 1)-BS_Pos(m, 1))^2+(RIS_Pos(n, 2)-BS_Pos(m, 2))^2+(RIS_y(end)-BS_y)^2); % distance between irs element shown by index and base station antenna as indicated by base station index
        g(n, m) = lambda/4/pi/Dis*exp(-1j*2*pi*Dis/lambda); % expression of channel for far field as frish formula (Type 1 channel)
    end
end

f = zeros(N, K+N*(R-1)); % channel between 2nd layer and user antenna , type 2 channel as given in [37] lemma 1
%% 2.Channel matrix betwenn user and IRS and IRS to IRS
for n = 1:N     % element index
    for k = 1:K % user antenna index
        f(n, k) = type2channel(lambda, RIS_Pos(n, 1), RIS_Pos(n, 2),User_Pos(k, 1), User_Pos(k, 2), RIS_y(1)-User_y);% channel betwen user and irs
    end
end
temp = zeros(N, N); % channel between layer 1 and layer 2
for n1 = 1:N     % element index of IRS1 layer
    for n2 = 1:N % element index of IRS2 layer
        temp(n1, n2) = type2channel(lambda, RIS_Pos(n1, 1), RIS_Pos(n1, 2),RIS_Pos(n2, 1), RIS_Pos(n2, 2), RIS_y(2)-RIS_y(1));% channel between IRS 1 layer and IRS layer 2
    end
end
f(:, K+1:end) = repmat(temp, 1, R-1); % making a single matrix of both the channel

%% Multi-Layer Transmit Beamformer Design or Proposed algorithm implementation
[SNR_MultiLayer_UC, UC_theta, UC_h, UC_Ratio1, UC_Ratio2, UC_SNR] = Proposed_Algorithm_multilayer(g, f, Sigma, Pmax, rpt, R, loss, N2, EUR_threshold);



%% Single-Layer Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Position

% also calculated above
%% 1.user generation
User_y = 0;             % y coordinate of user is at origin
User_Pos = zeros(K, 2); % (x,y) cordinate of K users
% generating (x,y) for all user antenna
for k = 1:K
    User_Pos(k, :) = [-(ceil(K/2)-1)*lambda/2-lambda/4+(k-1)*lambda/2, 0];% same as calculated above 
end

%% 2.single layer ris generation
RIS_y = Dis_RIS2User+Dis_Layer*(0:R-1); % y coordinate of RIS layer
RIS_Pos = zeros(N1(1)*N1(2), 2); % (x,y) coordinate of RIS layer
sideWidth = lambda/RISWidth; % diastance between layer
% generating position of each elements in multi  layer
for n = 1:N1(1)*N1(2)% element index of single layer case
    y = ceil(n/N1(2));
    x = n-(y-1)*N1(2);
    RIS_Pos(n, :) = [-(N1(2)/2-1)*sideWidth-sideWidth/2+(x-1)*sideWidth,(N1(1)/2-1)*sideWidth+sideWidth/2-(y-1)*sideWidth];
end


%% 3.base station antenna generation
BS_y = Dis_RIS2User*R+Dis_BS2RIS; % y coordinate of base station antennas (is same)
BS_Pos = zeros(M, 2);             % (x,y) coordinate of base station antennas
for k = 1:M
    BS_Pos(k, :) = [-(ceil(M/2)-1)*lambda/2-lambda/4+(k-1)*lambda/2, 0];
end

%% Generate Channel

N = N1(1)*N1(2);% total number of elements in layer

%% 1.channel between user and IRS
f = zeros(N, K+N*(R-1)); % channel between  layer and user antenna
for n = 1:N              % element index
    for k = 1:K          % user index
        f(n, k) = type2channel(lambda, RIS_Pos(n, 1), RIS_Pos(n, 2),User_Pos(k, 1), User_Pos(k, 2), RIS_y(1)-User_y);
    end
end
temp = zeros(N, N); % channel between layer 1 and layer 2
for n1 = 1:N 
    for n2 = 1:N
        temp(n1, n2) = type2channel(lambda, RIS_Pos(n1, 1), RIS_Pos(n1, 2),RIS_Pos(n2, 1), RIS_Pos(n2, 2), RIS_y(2)-RIS_y(1));
    end
end
f(:, K+1:end) = repmat(temp, 1, R-1);

%% 2.channel between irs layer and base station
g1 = zeros(N, M);
for n = 1:N     % element index
    for m = 1:M % base station index
        Dis = sqrt((RIS_Pos(n, 1)-BS_Pos(m, 1))^2+(RIS_Pos(n, 2)-BS_Pos(m, 2) )^2+(RIS_y(1)-BS_y)^2);%distance betwen indexed element and base station
        g1(n, m) = lambda/4/pi/Dis*exp(-1j*2*pi*Dis/lambda);% frish formula for far field
    end
end

%% Single-Layer beam forming and proposed algorithm implemented
[SNR_SingleLayer_CC, SNR_SingleLayer_UC, Tran_h, Tran_theta, Ref_SNR, Tran_SNR, Tran_Ratio] = Proposed_Algorithm_singlelayer(g1, f(:, 1:K), Sigma, Pmax, rpt, loss, N1, EUR_threshold);




%% No RIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Channel 
%% 1.channel between user and base station
g = zeros(K, M);
for k = 1:K
    for m = 1:M
        Dis =  sqrt((User_Pos(k, 1)-BS_Pos(m, 1))^2+(User_Pos(k, 2)-BS_Pos(m, 2))^2+(User_y-BS_y)^2);%distance between indexed element and base station antenna indexed
        g(k, m) = lambda/4/pi/Dis*exp(-1j*2*pi*Dis/lambda);% frish formula for far field channel
    end
end

%% beam forming
[v, ~] = eigs(g'*g, 1);% largest eigen value of g'g
a = g*v; % given in paper
w = a/norm(a); % given in paper
SNR_NoRIS = abs(v'*g'*w)^2/norm(v)^2/Sigma*Pmax; %given in paper

%% dB to magnitude conversion
SNR_SingleLayer_CC = pow2db(SNR_SingleLayer_CC);
SNR_MultiLayer_UC = pow2db(SNR_MultiLayer_UC);
SNR_SingleLayer_UC = pow2db(SNR_SingleLayer_UC);
SNR_NoRIS = pow2db(SNR_NoRIS);


%% Figure Plots%%%%%%%

%% Fig. 8. Detection SNR versus the maximum transmit power Pmax.    
figure;
hold on;
plot(Pmax_dB, SNR_MultiLayer_UC, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 9);
plot(Pmax_dB, SNR_SingleLayer_CC, 'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 9);
plot(Pmax_dB, SNR_SingleLayer_UC, 'LineWidth', 2, 'Marker', '^', 'MarkerSize', 9);
plot(Pmax_dB, SNR_NoRIS, 'LineWidth', 2, 'Marker', 's', 'MarkerSize', 9);
box on;
grid on;
title('Fig. 4. Detection SNR vs the maximum transmit power Pmax')
xlabel('${P}_{\rm{max}}$ (dB)', 'FontSize', 15, 'Interpreter', 'latex');
ylabel('SNR (dB)', 'FontSize', 15, 'Interpreter', 'latex');
legend('Multi-layer UC-RIS', 'Single-layer CC-RIS','Single-layer UC-RIS', 'No RIS', 'FontSize', 15, 'Location', 'NorthWest');

%% pattern
%% Fig. 9. Radiation pattern of different RISs.
N = size(UC_theta, 1); % number of elements in Transmit phase shift matrix of each layer of multi-layer
N2 = size(Tran_theta, 1); % number of elements in Transmit phase shift matrix of single layer case
idx = (0 : N-1)'; % indexing for multilayer
idx2 = (0 : N2-1)'; % indexing for single layer
AngleSet = linspace(0, 2*pi, 1E5); % set of angles that can be provided by each element with a resolutionof 10^-5
Hset = exp(1j * pi * idx * cos(AngleSet)); % exponential phase set for multilayer TPS from angle set
Hset2 = exp(1j * pi * idx2 * cos(AngleSet)); % exponential phase set for Single layer TPS from angle set
pat1 = diag(UC_theta(:, 1))*reshape(UC_h(:, :, 1), [], 1); % tps X column vector of 1st layer h
pat2 = diag(UC_theta(:, 2))*reshape(UC_h(:, :, 2), [], 1); % tps X column vector of 2st layer h
pat_Tran = diag(Tran_theta)*reshape(Tran_h, [], 1);

pat1 = reshape(pat1, [], 1);
pat2 = reshape(pat2, [], 1);
pat_Tran = reshape(pat_Tran, [], 1);
pat1_Transpose = reshape(pat1.', [], 1);
pat2_Transpose = reshape(pat2.', [], 1);
pat_Tran_Transpose = reshape(pat_Tran.', [], 1);

amp1 = abs(Hset'*pat1);
amp2 = abs(Hset'*pat2);
amp_Tran = abs(Hset2'*pat_Tran);
amp1_Transpose = abs(Hset'*pat1_Transpose);
amp2_Transpose = abs(Hset'*pat2_Transpose);
amp_Tran_Transpose = abs(Hset2'*pat_Tran_Transpose);

figure;
hold on;
box on;
plot(AngleSet/pi*180, amp2, 'Linewidth', 2);
plot(AngleSet/pi*180, amp_Tran, ':', 'Linewidth', 2);
plot(AngleSet/pi*180, amp1, '--', 'Linewidth', 2);
legend(['2^{nd} layer of multi-layer US-RIS'], 'Single-layer US-RIS', ['1^{st} layer of multi-layer US-RIS'], 'Fontsize', 15, 'Fontname', 'Times');
title('Fig. 5. Radiation pattern of different RISs.');
xlim([60 120]);
xticks(60:10:120);
xlabel('$\theta$ (deg)', 'Fontsize', 15, 'Interpreter', 'latex');
ylabel('Amplitude', 'Fontsize', 15, 'Interpreter', 'latex');
grid on;

%% convergence
%% Fig. 10. Detection SNR at Pmax = 0dBW against the number of iterations.
rpt = length(UC_SNR);
UC_SNR = 10*log10(UC_SNR);
Ref_SNR = 10*log10(Ref_SNR);
Tran_SNR = 10*log10(Tran_SNR);
figure;
hold on;
plot(0:rpt-1, UC_SNR, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 9);
plot(0:rpt-1, Ref_SNR, 'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 9);
plot(0:rpt-1, Tran_SNR, 'LineWidth', 2, 'Marker', '^', 'MarkerSize', 9);
box on;
grid on; 
xlim([0 rpt-1]);
title('Fig. 6. Detection SNR at Pmax = 0dBW against the number of iterations.')
xlabel('Iterations', 'FontSize', 15, 'Interpreter', 'latex');
ylabel('SNR (dB)', 'FontSize', 15, 'Interpreter', 'latex');
legend('Multi-layer UC-RIS', 'Single-layer CC-RIS', 'Single-layer UC-RIS', 'FontSize', 15, 'Location', 'SouthEast');

% Display EAR
[UC_Ratio1, UC_Ratio2, (UC_Ratio1+UC_Ratio2)/2, Tran_Ratio]

