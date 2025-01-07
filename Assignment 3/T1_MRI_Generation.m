% Author: Marwan Aridi
% Assignment 3
% This is a T1 Inversion file for the second part of the assignment
% Please see README file on how to run this code


% Load patient data and extract slice 90 for both patients
patient1_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient1.nii');
patient2_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient2.nii');
slice_90_patient1 = patient1_data(:, :, 90);
slice_90_patient2 = patient2_data(:, :, 90);

% Define TR and TI values for T1 Inversion Recovery table
TR_values = [1000, 2000];  % TR values in ms
TI_values = [50, 100, 250, 500, 750];  % TI values in ms

% Define T1, T2, and P value ranges for each tissue type
T1_GM_range = [1.62 - 0.26, 1.62 + 0.26]; T2_GM_range = [85 - 12, 85 + 12]; P_GM = 105.00;
T1_WM_range = [0.92 - 0.08, 0.92 + 0.08]; T2_WM_range = [80.5 - 9.75, 80.5 + 9.75]; P_WM = 80.00;
T1_CSF_range = [10.4 - 4.8, 10.4 + 4.8]; T2_CSF_range = [1055 - 472.5, 1055 + 472.5]; P_CSF = 150.00;

% Function to generate and display T1 Inversion Recovery table with noise
function generate_t1_ir_table(slice_data, TR_values, TI_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, patient_title)
    figure;
    tiledlayout(length(TI_values), length(TR_values), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(patient_title);
    
    for i = 1:length(TI_values)
        for j = 1:length(TR_values)
            TR = TR_values(j);
            TI = TI_values(i);

            % Generate 100 images with random T1 and T2 values and Gaussian noise
            for k = 1:100
                T1_GM = T1_GM_range(1) + (T1_GM_range(2) - T1_GM_range(1)) * rand;
                T2_GM = T2_GM_range(1) + (T2_GM_range(2) - T2_GM_range(1)) * rand;
                T1_WM = T1_WM_range(1) + (T1_WM_range(2) - T1_WM_range(1)) * rand;
                T2_WM = T2_WM_range(1) + (T2_WM_range(2) - T2_WM_range(1)) * rand;
                T1_CSF = T1_CSF_range(1) + (T1_CSF_range(2) - T1_CSF_range(1)) * rand;
                T2_CSF = T2_CSF_range(1) + (T2_CSF_range(2) - T2_CSF_range(1)) * rand;

                % Initialize signal intensity map
                SI_map = zeros(size(slice_data));
                
                % Create tissue masks based on intensity values in the slice
                GM_mask = (slice_data == 1);  WM_mask = (slice_data == 2);  CSF_mask = (slice_data == 3);

                % Calculate Signal Intensity for each tissue type
                SI_map(GM_mask) = P_GM * (1 - 2 * exp(-TI / (T1_GM * 1000)) + exp(-TR / (T1_GM * 1000)));
                SI_map(WM_mask) = P_WM * (1 - 2 * exp(-TI / (T1_WM * 1000)) + exp(-TR / (T1_WM * 1000)));
                SI_map(CSF_mask) = P_CSF * (1 - 2 * exp(-TI / (T1_CSF * 1000)) + exp(-TR / (T1_CSF * 1000)));

                % Add Gaussian noise with standard deviation = 5% of peak signal
                noise = 0.05 * max(SI_map(:)) * randn(size(SI_map));
                SI_map = SI_map + noise;
                
                % Store the first generated image
                if k == 1
                    selected_image = SI_map;
                end
            end
            
            % Display the selected image in the corresponding tile
            nexttile((i-1) * length(TR_values) + j);
            imshow(mat2gray(selected_image));  % Normalize and display
            title(sprintf('TR = %d, TI = %d', TR, TI), 'FontSize', 8);
        end
    end
end

% Generate and display T1 Inversion Recovery tables for Patient 1 and Patient 2
generate_t1_ir_table(slice_90_patient1, TR_values, TI_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, 'Patient 1 - T1 Inversion Recovery Table with Noise');
generate_t1_ir_table(slice_90_patient2, TR_values, TI_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF, 'Patient 2 - T1 Inversion Recovery Table with Noise');

% Function to generate T1 Inversion Recovery values
function T1_IR_values = generate_t1_ir_table_values(slice_data, TR_values, TI_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF)
    T1_IR_values = zeros(length(TI_values), length(TR_values), 3); % 3 tissues: GM, WM, CSF

    for i = 1:length(TI_values)
        for j = 1:length(TR_values)
            TR = TR_values(j);
            TI = TI_values(i);

            % Sample T1 and T2 values for each tissue
            T1_GM = T1_GM_range(1) + (T1_GM_range(2) - T1_GM_range(1)) * rand;
            T2_GM = T2_GM_range(1) + (T2_GM_range(2) - T2_GM_range(1)) * rand;
            T1_WM = T1_WM_range(1) + (T1_WM_range(2) - T1_WM_range(1)) * rand;
            T2_WM = T2_WM_range(1) + (T2_WM_range(2) - T2_WM_range(1)) * rand;
            T1_CSF = T1_CSF_range(1) + (T1_CSF_range(2) - T1_CSF_range(1)) * rand;
            T2_CSF = T2_CSF_range(1) + (T2_CSF_range(2) - T2_CSF_range(1)) * rand;

            % Calculate Signal Intensity for each tissue type
            T1_IR_values(i, j, 1) = P_GM * (1 - 2 * exp(-TI / (T1_GM * 1000)) + exp(-TR / (T1_GM * 1000)));
            T1_IR_values(i, j, 2) = P_WM * (1 - 2 * exp(-TI / (T1_WM * 1000)) + exp(-TR / (T1_WM * 1000)));
            T1_IR_values(i, j, 3) = P_CSF * (1 - 2 * exp(-TI / (T1_CSF * 1000)) + exp(-TR / (T1_CSF * 1000)));
        end
    end
end

% Retrieve and display T1 Inversion Recovery values for Patient 1
T1_IR_values_patient1 = generate_t1_ir_table_values(slice_90_patient1, TR_values, TI_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF);
disp('T1 Inversion Recovery Values for Patient 1:');
display_t1_ir_values(T1_IR_values_patient1, TR_values, TI_values);

% Retrieve and display T1 Inversion Recovery values for Patient 2
T1_IR_values_patient2 = generate_t1_ir_table_values(slice_90_patient2, TR_values, TI_values, T1_GM_range, T2_GM_range, T1_WM_range, T2_WM_range, T1_CSF_range, T2_CSF_range, P_GM, P_WM, P_CSF);
disp('T1 Inversion Recovery Values for Patient 2:');
display_t1_ir_values(T1_IR_values_patient2, TR_values, TI_values);

% Helper function to display T1 Inversion Recovery values
function display_t1_ir_values(T1_IR_values, TR_values, TI_values)
    tissue_names = {'GM', 'WM', 'CSF'};
    for tissue = 1:3
        disp(['Tissue: ', tissue_names{tissue}]);
        tissue_data = squeeze(T1_IR_values(:, :, tissue));
        TR_column_names = arrayfun(@(tr) sprintf('TR_%d', tr), TR_values, 'UniformOutput', false);
        disp(array2table(tissue_data, 'VariableNames', TR_column_names, 'RowNames', arrayfun(@(ti) sprintf('TI_%d', TI_values(ti)), 1:length(TI_values), 'UniformOutput', false)));
    end
end
