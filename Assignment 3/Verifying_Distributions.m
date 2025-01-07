% Author: Marwan Aridi
% Assignment 3
% This is the code for the last part of this assignment
% Please see README file on how to run this code


% Define ranges for T1, T2, and P values for each tissue type (from Table 1)
% Gray Matter (GM)
T1_GM_range = [1.62 - 0.26, 1.62 + 0.26]; 
T2_GM_range = [85 - 12, 85 + 12];
P_GM = 105.00;

% White Matter (WM)
T1_WM_range = [0.92 - 0.08, 0.92 + 0.08];
T2_WM_range = [80.5 - 9.75, 80.5 + 9.75];
P_WM = 80.00;

% Cerebrospinal Fluid (CSF)
T1_CSF_range = [10.4 - 4.8, 10.4 + 4.8];
T2_CSF_range = [1055 - 472.5, 1055 + 472.5];
P_CSF = 150.00;

% Number of images per tissue type for distribution analysis
num_images = 100;

% Initialize structures to store T1, T2, P, and SNR values for each tissue
tissue_types = {'GM', 'WM', 'CSF'};
T1_values = struct('GM', [], 'WM', [], 'CSF', []);
T2_values = struct('GM', [], 'WM', [], 'CSF', []);
P_values = struct('GM', [], 'WM', [], 'CSF', []);
SNR_values = struct('GM', [], 'WM', [], 'CSF', []);

% Standard deviation for Gaussian noise as a percentage of peak signal
noise_std_ratio = 0.05;

% Loop to generate values for 100 synthetic images per tissue type
for k = 1:num_images
    % Randomly sample T1 and T2 values within specified ranges for each tissue
    T1_GM = T1_GM_range(1) + (T1_GM_range(2) - T1_GM_range(1)) * rand;
    T2_GM = T2_GM_range(1) + (T2_GM_range(2) - T2_GM_range(1)) * rand;
    T1_WM = T1_WM_range(1) + (T1_WM_range(2) - T1_WM_range(1)) * rand;
    T2_WM = T2_WM_range(1) + (T2_WM_range(2) - T2_WM_range(1)) * rand;
    T1_CSF = T1_CSF_range(1) + (T1_CSF_range(2) - T1_CSF_range(1)) * rand;
    T2_CSF = T2_CSF_range(1) + (T2_CSF_range(2) - T2_CSF_range(1)) * rand;

    % Store T1 and T2 values for each tissue type
    T1_values.GM = [T1_values.GM, T1_GM];
    T2_values.GM = [T2_values.GM, T2_GM];
    T1_values.WM = [T1_values.WM, T1_WM];
    T2_values.WM = [T2_values.WM, T2_WM];
    T1_values.CSF = [T1_values.CSF, T1_CSF];
    T2_values.CSF = [T2_values.CSF, T2_CSF];

    % Store P values for each tissue type (constant for each tissue)
    P_values.GM = [P_values.GM, P_GM];
    P_values.WM = [P_values.WM, P_WM];
    P_values.CSF = [P_values.CSF, P_CSF];

    % Calculate noisy signal for each tissue and SNR
    % Adding Gaussian noise proportional to P value
    noise_GM = noise_std_ratio * P_GM * randn;
    noise_WM = noise_std_ratio * P_WM * randn;
    noise_CSF = noise_std_ratio * P_CSF * randn;
    
    % Signal with noise for each tissue
    noisy_signal_GM = P_GM + noise_GM;
    noisy_signal_WM = P_WM + noise_WM;
    noisy_signal_CSF = P_CSF + noise_CSF;

    % Calculate SNR for each tissue as mean signal / noise standard deviation
    SNR_GM = noisy_signal_GM / abs(noise_GM);
    SNR_WM = noisy_signal_WM / abs(noise_WM);
    SNR_CSF = noisy_signal_CSF / abs(noise_CSF);

    % Store SNR values for each tissue type
    SNR_values.GM = [SNR_values.GM, SNR_GM];
    SNR_values.WM = [SNR_values.WM, SNR_WM];
    SNR_values.CSF = [SNR_values.CSF, SNR_CSF];
end

% Calculate mean and standard deviation for each parameter for each tissue
results = struct();
for i = 1:length(tissue_types)
    tissue = tissue_types{i};
    results.(tissue).T1_mean = mean(T1_values.(tissue));
    results.(tissue).T1_std = std(T1_values.(tissue));
    results.(tissue).T2_mean = mean(T2_values.(tissue));
    results.(tissue).T2_std = std(T2_values.(tissue));
    results.(tissue).P_mean = mean(P_values.(tissue));
    results.(tissue).P_std = std(P_values.(tissue));
    results.(tissue).SNR_mean = mean(SNR_values.(tissue));
    results.(tissue).SNR_std = std(SNR_values.(tissue));
end

% Display results in a detailed format
disp('Results of T1, T2, P, and SNR Analysis for Gray Matter (GM):');
disp(results.GM);

disp('Results of T1, T2, P, and SNR Analysis for White Matter (WM):');
disp(results.WM);

disp('Results of T1, T2, P, and SNR Analysis for Cerebrospinal Fluid (CSF):');
disp(results.CSF);
