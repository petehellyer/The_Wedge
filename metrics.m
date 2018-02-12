function [phi] = metrics(theta)
% Calculate synchronisation index (Sync), metastability (Meta), and
% instantaneous synchrony (phi), for phase data theta
N = size(theta,2);
Tmax = size(theta,1);
% Calculate synchrony
phi = zeros(Tmax,1);
for j = 1:N
    x = theta(:,j);
    phi = phi + exp(x*sqrt(-1));
end
phi = phi/N;
end