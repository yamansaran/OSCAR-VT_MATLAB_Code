% Sigma values (x-axis)
sigma = 1:10;

% -------------------------
% NO SNS DATA (meters)
% -------------------------
R_no = [101 202 302 403 504 605 705 806 907 1010];
T_no = [948 1900 2840 3790 4740 5690 6630 7580 8530 9480];
N_no = [72 144 215 287 359 431 503 574 646 718];

% Compute RSS (No SNS)
RSS_no = sqrt(R_no.^2 + T_no.^2 + N_no.^2);

% -------------------------
% WITH SNS DATA (meters)
% -------------------------
R_sns = [3 6 9 12 15 18 21 24 27 30];
T_sns = [3 6 9 12 15 18 21 24 27 30];
N_sns = [3 6 9 12 15 18 21 24 27 30];

% Compute RSS (SNS)
RSS_sns = sqrt(R_sns.^2 + T_sns.^2 + N_sns.^2);

% =========================
% FIGURE 1: NO SNS
% =========================
figure;
plot(sigma, R_no, '-o', 'LineWidth', 1.5); hold on;
plot(sigma, T_no, '-s', 'LineWidth', 1.5);
plot(sigma, N_no, '-^', 'LineWidth', 1.5);
plot(sigma, RSS_no, '-d', 'LineWidth', 1.5);

grid on;
xlabel('\sigma Level');
ylabel('Distance (m)');
title('Position Uncertainty vs Sigma (No Sensor Fusion)');
legend('Radial (R)', 'Along-Track (T)', 'Cross-Track (N)', 'RSS', 'Location', 'northwest');

% =========================
% FIGURE 2: WITH SNS
% =========================
figure;
plot(sigma, R_sns, '-o', 'LineWidth', 1.5); hold on;
plot(sigma, T_sns, '-s', 'LineWidth', 1.5);
plot(sigma, N_sns, '-^', 'LineWidth', 1.5);
plot(sigma, RSS_sns, '-d', 'LineWidth', 1.5);

grid on;
xlabel('\sigma Level');
ylabel('Distance (m)');
title('Position Uncertainty vs Sigma (With Sensor Fusion)');
legend('Radial (R)', 'Along-Track (T)', 'Cross-Track (N)', 'RSS', 'Location', 'northwest');