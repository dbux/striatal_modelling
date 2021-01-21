function[E] = gen_e
% Generates the value E for each connection type as specified on page 9 of
% Humphries, Wood & Gurney (2010)

% Set maximum distance to greater than striatal edge length to allow for 
% corner-to-corner checks
d = 1:2000;

% Values from page 10, table 4
alpha = [0.511, -0.921, -0.695, 1.322];
beta  = [1.033, 1.033, 1.38, 2.4];
gamma = [0.042, 0.042, 0.057, 0.016];
delta = [26.8, 26.8, 15.6, 43.3];
eta   = [0.0039, 0.0039, 0.0036, 0.0029];

% Page 9, equation 20
for c = 1:numel(alpha) 
    E(c,:) = exp(-alpha(c) - beta(c) .* (1 - exp(-gamma(c) .* (d - delta(c)))) .* exp(eta(c) .* d));
end

% % Optionally generate plot of E(c)
% % Page 9, figure 7
% figure(1); clf;
% semilogy(E(1,:), 'r', 'Linewidth', 2)
% xlabel('Distance d_s between somas (µm)')
% ylabel('Expected number of contacts')
% hold on
% grid on
% axis([0 600 0.0001 10])
% semilogy(E(2,:), 'b', 'Linewidth', 2)
% semilogy(E(3,:), 'g', 'Linewidth', 2)
% semilogy(E(4,:), 'k', 'Linewidth', 2)
% legend('MSN-MSN', 'FSI-MSN', 'FSI-FSI', 'FSI gap')