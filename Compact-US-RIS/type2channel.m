%% ECN 620 TERM PAPER CODE
%% Nitish 21531009, Manmohan 21531006 

%Power Scaling Laws and Near-Field Behaviors of Massive MIMO and Intelligent Reflecting Surfaces Emil Björnson, Senior Member, IEEE, Luca Sanguinetti, Senior Member, IEEE
% As given in [37], lemma 1 of above paper
% same gain expression
function h = type2channel(lambda, xn, zn, xt, zt, d) % variable as given in 
X = [lambda/4+xn-xt, lambda/4-xn+xt];
Z = [lambda/4+zn-zt, lambda/4-zn+zt];
h = 0;
for i = 1:2
    for j = 1:2
        x = X(i);
        z = Z(j);
        l = sqrt(x^2+z^2+d^2);
        h = h+x*z*d/3/(z^2+d^2)/l+2/3*atan(x*z/d/l);
    end
end
h = h/4/pi;
h = sqrt(h);% square root of gain
h = h*exp(-1j*2*pi*d/lambda); %friss formula for above gain
end