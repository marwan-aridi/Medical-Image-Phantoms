% Author: Marwan Aridi
% Assignment 3
% This is a GRE file for the second part of the assignment
% Please see README file on how to run this code


% Load patient data and extract slice 90 for both patients
patient1_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient1.nii');
patient2_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient2.nii');
slice_90_patient1 = patient1_data(:, :, 90);
slice_90_patient2 = patient2_data(:, :, 90);

% Define TR and flip angle (α) values for GRE table
TR_values = [25, 50, 100, 200];           % TR values in ms
alpha_values = [15, 30, 45, 60, 90];      % Flip angle values in degrees
TE_fixed = 5;                             % Fixed TE for GRE in ms

% Define T1, T2, and P value ranges for each tissue type
T1_GM_range = [1.62 - 0.26, 1.62 + 0.26]; T2_GM_range = [85 - 12, 85 + 12]; P_GM = 105.00;
T1_WM_range = [0.92 - 0.08, 0.92 + 0.08]; T2_WM_range = [80.5 - 9.75, 80.5 + 9.75]; P_WM = 80.00;
T1_CSF_range = [10.4 - 4.8, 10.4 + 4.8]; T2_CSF_range = [1055 - 472.5, 1055 + 472.5]; P_CSF = 150.00;

% Function to generate and display GRE table with noise
function generate_gre_table(slice_data, TR_values, alpha_values, TE_fixed, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, patient_title)
    figure;
    tiledlayout(length(alpha_values), length(TR_values), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(patient_title);
    
    for i = 1:length(alpha_values)
        for j = 1:length(TR_values)
            TR = TR_values(j);
            alpha = alpha_values(i);

            % Generate 100 images with random T1 and T2 values and Gaussian noise
            for k = 1:100
                T1_GM = T1_GM_range(1) + (T1_GM_range(2) - T1_GM_range(1)) * rand;
                T2_GM = T2_GM_range(1) + (T2_GM_range(2) - T2_GM_range(1)) * rand;
                T1_WM = T1_WM_range(1) + (T1_WM_range(2) - T1_WM_range(1)) * rand;
                T2_WM = T2_WM_range(1) + (T2_WM_range(2) - T2_WM_range(1)) * rand;
                T1_CSF = T1_CSF_range(1) + (T1_CSF_range(2) - T1_CSF_range(1)) * rand;
                T2_CSF = T2_CSF_range(1) + (T2_CSF_range(2) - T2_CSF_range(1)) * rand;

                SI_map = zeros(size(slice_data));
                GM_mask = (slice_data == 1); WM_mask = (slice_data == 2); CSF_mask = (slice_data == 3);
                alpha_rad = deg2rad(alpha);

                % Calculate GRE Signal Intensity for each tissue
                SI_map(GM_mask) = P_GM * (1 - exp(-TR / (T1_GM * 1000))) * sin(alpha_rad) * exp(-TE_fixed / T2_GM) / ...
                                  (1 - cos(alpha_rad) * exp(-TR / (T1_GM * 1000)));
                SI_map(WM_mask) = P_WM * (1 - exp(-TR / (T1_WM * 1000))) * sin(alpha_rad) * exp(-TE_fixed / T2_WM) / ...
                                  (1 - cos(alpha_rad) * exp(-TR / (T1_WM * 1000)));
                SI_map(CSF_mask) = P_CSF * (1 - exp(-TR / (T1_CSF * 1000))) * sin(alpha_rad) * exp(-TE_fixed / T2_CSF) / ...
                                   (1 - cos(alpha_rad) * exp(-TR / (T1_CSF * 1000)));
                
                % Add Gaussian noise with standard deviation = 5% of peak signal
                noise = 0.05 * max(SI_map(:)) * randn(size(SI_map));
                SI_map = SI_map + noise;

                % Store this version of SI_map if needed for analysis
                if k == 1
                    selected_image = SI_map;
                end
            end
            
            % Display the selected image in the corresponding tile
            nexttile((i-1) * length(TR_values) + j);
            imshow(mat2gray(selected_image));  % Normalize and display
            title(sprintf('TR = %d, α = %d°', TR, alpha), 'FontSize', 8);
        end
    end
end

% Generate GRE tables for both patients
generate_gre_table(slice_90_patient1, TR_values, alpha_values, TE_fixed, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, 'Patient 1 - GRE Table with Noise');
generate_gre_table(slice_90_patient2, TR_values, alpha_values, TE_fixed, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, 'Patient 2 - GRE Table with Noise');

% Function to generate GRE values for each tissue type
function GRE_values = generate_gre_values(slice_data, TR_values, alpha_values, TE_fixed, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF)
    GRE_values = zeros(length(alpha_values), length(TR_values), 3); % 3 tissues: GM, WM, CSF

    for i = 1:length(alpha_values)
        for j = 1:length(TR_values)
            TR = TR_values(j);
            alpha = alpha_values(i);
            alpha_rad = deg2rad(alpha);

            % Sample T1 and T2 values within specified ranges for each tissue
            T1_GM = T1_GM_range(1) + (T1_GM_range(2) - T1_GM_range(1)) * rand;
            T2_GM = T2_GM_range(1) + (T2_GM_range(2) - T2_GM_range(1)) * rand;
            T1_WM = T1_WM_range(1) + (T1_WM_range(2) - T1_WM_range(1)) * rand;
            T2_WM = T2_WM_range(1) + (T2_WM_range(2) - T2_WM_range(1)) * rand;
            T1_CSF = T1_CSF_range(1) + (T1_CSF_range(2) - T1_CSF_range(1)) * rand;
            T2_CSF = T2_CSF_range(1) + (T2_CSF_range(2) - T2_CSF_range(1)) * rand;

            % Calculate Signal Intensity using GRE formula
            GRE_values(i, j, 1) = P_GM * (1 - exp(-TR / (T1_GM * 1000))) * sin(alpha_rad) * exp(-TE_fixed / T2_GM) / ...
                                  (1 - cos(alpha_rad) * exp(-TR / (T1_GM * 1000)));
            GRE_values(i, j, 2) = P_WM * (1 - exp(-TR / (T1_WM * 1000))) * sin(alpha_rad) * exp(-TE_fixed / T2_WM) / ...
                                  (1 - cos(alpha_rad) * exp(-TR / (T1_WM * 1000)));
            GRE_values(i, j, 3) = P_CSF * (1 - exp(-TR / (T1_CSF * 1000))) * sin(alpha_rad) * exp(-TE_fixed / T2_CSF) / ...
                                  (1 - cos(alpha_rad) * exp(-TR / (T1_CSF * 1000)));
        end
    end
end

% Retrieve and display GRE values for both patients in LaTeX format
GRE_values_patient1 = generate_gre_values(slice_90_patient1, TR_values, alpha_values, TE_fixed, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF);
GRE_values_patient2 = generate_gre_values(slice_90_patient2, TR_values, alpha_values, TE_fixed, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF);

disp('GRE Values for Patient 1:');
display_gre_values_latex(GRE_values_patient1, TR_values, alpha_values);

disp('GRE Values for Patient 2:');
display_gre_values_latex(GRE_values_patient2, TR_values, alpha_values);

% Helper function to display GRE values in LaTeX format
function display_gre_values_latex(GRE_values, TR_values, alpha_values)
    tissue_names = {'GM', 'WM', 'CSF'};
    for tissue = 1:3
        disp(['\textbf{', tissue_names{tissue}, '}']);
        for i = 1:length(alpha_values)
            fprintf('%d & ', alpha_values(i));
            for j = 1:length(TR_values)
                fprintf('%.4f ', GRE_values(i, j, tissue));
                if j < length(TR_values)
                    fprintf('& ');
                else
                    fprintf('\\\\ \n');
                end
            end
        end
    end
end
