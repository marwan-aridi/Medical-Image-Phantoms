% Author: Marwan Aridi
% Assignment 3
% This is a Spin Echo file for the second part of the assignment
% Please see README file on how to run this code


% Load patient data and extract slice 90 for both patients
patient1_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient1.nii');
patient2_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient2.nii');
slice_90_patient1 = patient1_data(:, :, 90);
slice_90_patient2 = patient2_data(:, :, 90);

% Define TR and TE values for Spin Echo table
TR_values = [250, 500, 750, 1000, 2000];   % TR values in ms
TE_values = [20, 40, 60, 80];              % TE values in ms

% Define T1, T2, and P value ranges for each tissue type (from Table 1)
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

% Function to generate and display Spin Echo table with noise for a given patient slice
function generate_spin_echo_table(slice_data, TR_values, TE_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, patient_title)
    figure;
    tiledlayout(length(TR_values), length(TE_values), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(patient_title);
    
    for i = 1:length(TR_values)
        for j = 1:length(TE_values)
            TR = TR_values(i);
            TE = TE_values(j);

            % Generate 100 images with random T1 and T2 values and Gaussian noise
            for k = 1:100
                % Randomly sample T1 and T2 within specified ranges for each tissue
                T1_GM = T1_GM_range(1) + (T1_GM_range(2) - T1_GM_range(1)) * rand;
                T2_GM = T2_GM_range(1) + (T2_GM_range(2) - T2_GM_range(1)) * rand;
                T1_WM = T1_WM_range(1) + (T1_WM_range(2) - T1_WM_range(1)) * rand;
                T2_WM = T2_WM_range(1) + (T2_WM_range(2) - T2_WM_range(1)) * rand;
                T1_CSF = T1_CSF_range(1) + (T1_CSF_range(2) - T1_CSF_range(1)) * rand;
                T2_CSF = T2_CSF_range(1) + (T2_CSF_range(2) - T2_CSF_range(1)) * rand;

                % Initialize the signal intensity map for the patient
                SI_map = zeros(size(slice_data));
                
                % Create tissue masks based on intensity values in the slice
                GM_mask = (slice_data == 1);  % Gray Matter
                WM_mask = (slice_data == 2);  % White Matter
                CSF_mask = (slice_data == 3); % Cerebrospinal Fluid

                % Calculate Signal Intensity for each tissue type and combine
                % Gray Matter
                SI_map(GM_mask) = P_GM * (1 - exp(-TR / (T1_GM * 1000))) * exp(-TE / T2_GM);
                
                % White Matter
                SI_map(WM_mask) = P_WM * (1 - exp(-TR / (T1_WM * 1000))) * exp(-TE / T2_WM);
                
                % Cerebrospinal Fluid
                SI_map(CSF_mask) = P_CSF * (1 - exp(-TR / (T1_CSF * 1000))) * exp(-TE / T2_CSF);
                
                % Add Gaussian noise with standard deviation = 5% of peak signal
                peak_signal = max(SI_map(:));
                noise = 0.05 * peak_signal * randn(size(SI_map));
                SI_map = SI_map + noise;
                
                % Store this version of SI_map if needed for analysis
                if k == 1
                    % Select the first generated image (or adjust as needed)
                    selected_image = SI_map;
                end
            end
            
            % Display the selected image in the corresponding tile
            nexttile((i-1) * length(TE_values) + j);
            imshow(mat2gray(selected_image));  % Normalize and display
            title(sprintf('TR = %d, TE = %d', TR, TE), 'FontSize', 8);
        end
    end
end

% Generate Spin Echo table for Patient 1
generate_spin_echo_table(slice_90_patient1, TR_values, TE_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, 'Patient 1 - Spin Echo Table with Noise');

% Generate Spin Echo table for Patient 2
generate_spin_echo_table(slice_90_patient2, TR_values, TE_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, 'Patient 2 - Spin Echo Table with Noise');

% Function to calculate mean signal intensity for each tissue type
function mean_intensities = calculate_mean_intensities(slice_data, TR_values, TE_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF)
    mean_intensities = struct('GM', [], 'WM', [], 'CSF', []);
    
    for i = 1:length(TR_values)
        for j = 1:length(TE_values)
            TR = TR_values(i);
            TE = TE_values(j);
            
            % Generate 100 images with random T1 and T2 values and Gaussian noise
            gm_signals = [];
            wm_signals = [];
            csf_signals = [];
            
            for k = 1:100
                % Randomly sample T1 and T2 within specified ranges for each tissue
                T1_GM = T1_GM_range(1) + (T1_GM_range(2) - T1_GM_range(1)) * rand;
                T2_GM = T2_GM_range(1) + (T2_GM_range(2) - T2_GM_range(1)) * rand;
                T1_WM = T1_WM_range(1) + (T1_WM_range(2) - T1_WM_range(1)) * rand;
                T2_WM = T2_WM_range(1) + (T2_WM_range(2) - T2_WM_range(1)) * rand;
                T1_CSF = T1_CSF_range(1) + (T1_CSF_range(2) - T1_CSF_range(1)) * rand;
                T2_CSF = T2_CSF_range(1) + (T2_CSF_range(2) - T2_CSF_range(1)) * rand;

                % Create tissue masks based on intensity values in the slice
                GM_mask = (slice_data == 1);  % Gray Matter
                WM_mask = (slice_data == 2);  % White Matter
                CSF_mask = (slice_data == 3); % Cerebrospinal Fluid

                % Calculate Signal Intensity for each tissue type
                gm_intensity = P_GM * (1 - exp(-TR / (T1_GM * 1000))) * exp(-TE / T2_GM);
                wm_intensity = P_WM * (1 - exp(-TR / (T1_WM * 1000))) * exp(-TE / T2_WM);
                csf_intensity = P_CSF * (1 - exp(-TR / (T1_CSF * 1000))) * exp(-TE / T2_CSF);

                % Add Gaussian noise with standard deviation = 5% of signal
                gm_signals = [gm_signals; gm_intensity * (1 + 0.05 * randn)];
                wm_signals = [wm_signals; wm_intensity * (1 + 0.05 * randn)];
                csf_signals = [csf_signals; csf_intensity * (1 + 0.05 * randn)];
            end
            
            % Compute the mean signal for each tissue and store
            mean_intensities.GM(i, j) = mean(gm_signals);
            mean_intensities.WM(i, j) = mean(wm_signals);
            mean_intensities.CSF(i, j) = mean(csf_signals);
        end
    end
end

% Calculate and display mean signal intensities for Patient 1
mean_intensities_patient1 = calculate_mean_intensities(slice_90_patient1, TR_values, TE_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF);
disp('Gray Matter (GM) Signal Intensity Values for Patient 1:');
disp(mean_intensities_patient1.GM);
disp('White Matter (WM) Signal Intensity Values for Patient 1:');
disp(mean_intensities_patient1.WM);
disp('Cerebrospinal Fluid (CSF) Signal Intensity Values for Patient 1:');
disp(mean_intensities_patient1.CSF);

% Calculate and display mean signal intensities for Patient 2
mean_intensities_patient2 = calculate_mean_intensities(slice_90_patient2, TR_values, TE_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF);
disp('Gray Matter (GM) Signal Intensity Values for Patient 2:');
disp(mean_intensities_patient2.GM);
disp('White Matter (WM) Signal Intensity Values for Patient 2:');
disp(mean_intensities_patient2.WM);
disp('Cerebrospinal Fluid (CSF) Signal Intensity Values for Patient 2:');
disp(mean_intensities_patient2.CSF);
