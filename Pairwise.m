% =========================================================
% Analytic Hierarchy Process (AHP) Interactive Comparison
% =========================================================
clear; clc;

% ---- Define your parameters ----
params = {'Mass (kg/m^2)', 'Power (W/m^2)', 'Thickness (mm)','EOL/PA5'};
n = length(params);clc


fprintf('\n=== AHP Interactive Pairwise Comparison ===\n');
fprintf('You will compare parameters using the Saaty 1–9 scale:\n');
fprintf(' 1 = Equal importance\n 3 = Moderate importance\n 5 = Strong importance\n 7 = Very strong\n 9 = Extreme\n');
fprintf(' Use fractions (e.g., 1/3) if the second is more important.\n\n');

% ---- Initialize matrix ----
A = ones(n);

% ---- Collect pairwise inputs ----
for i = 1:n-1
    for j = i+1:n
        prompt = sprintf('How much more important is "%s" than "%s"? ', params{i}, params{j});
        val = input(prompt);
        A(i,j) = val;
        A(j,i) = 1/val;
    end
end

% ---- Eigenvector method ----
[V, D] = eig(A);
[lambda_max, idx] = max(real(diag(D)));
w = real(V(:, idx));
w = w / sum(w);

% ---- Geometric mean method ----
geom_mean = prod(A, 2).^(1/n);
w_geom = geom_mean / sum(geom_mean);

% ---- Consistency Check ----
CI = (lambda_max - n) / (n - 1);
RI_table = [0 0 0.58 0.90 1.12 1.24 1.32 1.41 1.45 1.49];
RI = RI_table(n);
CR = CI / RI;

% ---- Display Results ----
fprintf('\n=== AHP Results ===\n');
fprintf('Lambda_max: %.4f\n', lambda_max);
fprintf('CI: %.4f\n', CI);
fprintf('RI: %.2f\n', RI);
fprintf('CR: %.4f\n', CR);

fprintf('\nParameter Weights (Eigenvector Method):\n');
for i = 1:n
    fprintf('  %-15s : %.4f\n', params{i}, w(i));
end

fprintf('\nParameter Weights (Geometric Mean Method):\n');
for i = 1:n
    fprintf('  %-15s : %.4f\n', params{i}, w_geom(i));
end

if CR < 0.10
    fprintf('\nConsistency Ratio (%.4f) is acceptable\n', CR);
else
    fprintf('\nConsistency Ratio (%.4f) is too high ️ — review comparisons.\n', CR);
end

% ---- Optional: display pairwise matrix ----
disp(' ');
disp('Pairwise Comparison Matrix:');
disp(array2table(A, 'VariableNames', params, 'RowNames', params));
