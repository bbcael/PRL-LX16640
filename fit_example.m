close all; clear all; clc; %clear workspace
load X % load data

X = X(~isnan(X)); X = X(X>0); % remove nonpositive values

% now fit log-skew-normal
ee = sort(log10(X)-min(log10(X)));
[y x] = ecdf(ee);
x = x(2:end); y = y(2:end);
q = unique(ee);
[a0b0c0] = [1.9319 18.4152 2.8145]; % initial parameter guesses
a0 = a0b0c0(1);
b0 = a0b0c0(2);
c0 = a0b0c0(3);
snr = 100; % choose signal to noise ratio for jumps in parameter space
lsncdf = normcdf((10.*log(q)-b0)/c0)-2.*tfn((10.*log(q)-b0)/c0,a0);
l0 = max(lsncdf-y);
u0 = max(y-lsncdf);
ul0 = u0+l0;
UL = [];
for i = 1:10000; % iteratively jump in parameter space to improve fit
    i
    a1 = awgn(a0,snr);
    b1 = awgn(b0,snr);
    c1 = awgn(c0,snr);
    lsncdf = normcdf((10.*log(q)-b1)/c1)-2.*tfn((10.*log(q)-b1)/c1,a1);
    l = max(lsncdf-y);
    u = max(y-lsncdf);
    if l+u <= ul0; % update if guess is better
        a0 = a1; b0 = b1; c0 = c1; ul0 = l+u; u0 = u; l0 = l;
    end
    UL(i) = u0+l0;
end
plot(1:length(UL),UL) % plot goodness-of-fit difference over iterations
