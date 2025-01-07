% Author: Marwan Aridi
% Assignment 3
% This is a GRE file for the first part of the assignment
% Please see README file on how to run this code

% Load patient data and extract slice 90
patient1_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient1.nii');
slice_90_patient1 = patient1_data(:, :, 90);

% Define TR and flip angle values for Gradient Recalled Echo table
TR_values = [25, 50, 100, 200];   % TR values in ms
alpha_values = [15, 30, 45, 60, 90];  % Flip angles in degrees
TE_GRE = 5;  % Fixed TE value for Gradient Recalled Echo

% Define T1, T2, and P values for each tissue type (GM, WM, CSF)
% Gray Matter (GM)
T1_GM = 1.62; T2_GM = 85.00; P_GM = 105.00;

% White Matter (WM)
T1_WM = 0.92; T2_WM = 80.50; P_WM = 80.00;

% Cerebrospinal Fluid (CSF)
T1_CSF = 10.4; T2_CSF = 1055.00; P_CSF = 150.00;

% Initialize arrays to store the calculated signal intensity values
SI_values_GM = zeros(length(alpha_values), length(TR_values));
SI_values_WM = zeros(length(alpha_values), length(TR_values));
SI_values_CSF = zeros(length(alpha_values), length(TR_values));

% Set up figure for displaying images
figure;
tiledlayout(length(alpha_values), length(TR_values), 'TileSpacing', 'compact', 'Padding', 'compact');
title('Patient 1 - Gradient Recalled Echo Table');

% Calculate the Signal Intensity for each combination of TR and alpha
for i = 1:length(alpha_values)
    for j = 1:length(TR_values)
        TR = TR_values(j);
        alpha = alpha_values(i);
        alpha_rad = deg2rad(alpha); % Convert flip angle to radians

        % Initialize the signal intensity map for Patient 1
        SI_map = zeros(size(slice_90_patient1));
        
        % Create tissue masks based on intensity values in the slice
        GM_mask = (slice_90_patient1 == 1);  % Gray Matter
        WM_mask = (slice_90_patient1 == 2);  % White Matter
        CSF_mask = (slice_90_patient1 == 3); % Cerebrospinal Fluid

        % Calculate Signal Intensity for each tissue type and save to arrays
        % Gray Matter
        SI_values_GM(i, j) = P_GM * (1 - exp(-TR / (T1_GM * 1000))) * sin(alpha_rad) * exp(-TE_GRE / T2_GM) / ...
                             (1 - cos(alpha_rad) * exp(-TR / (T1_GM * 1000)));
        SI_map(GM_mask) = SI_values_GM(i, j);
        
        % White Matter
        SI_values_WM(i, j) = P_WM * (1 - exp(-TR / (T1_WM * 1000))) * sin(alpha_rad) * exp(-TE_GRE / T2_WM) / ...
                             (1 - cos(alpha_rad) * exp(-TR / (T1_WM * 1000)));
        SI_map(WM_mask) = SI_values_WM(i, j);
        
        % Cerebrospinal Fluid
        SI_values_CSF(i, j) = P_CSF * (1 - exp(-TR / (T1_CSF * 1000))) * sin(alpha_rad) * exp(-TE_GRE / T2_CSF) / ...
                              (1 - cos(alpha_rad) * exp(-TR / (T1_CSF * 1000)));
        SI_map(CSF_mask) = SI_values_CSF(i, j);
        
        % Display the image in the corresponding tile
        nexttile((i-1) * length(TR_values) + j);
        imshow(mat2gray(SI_map));  % Normalize and display
        title(sprintf('TR = %d, α = %d°', TR, alpha), 'FontSize', 8);
    end
end

% Display the calculated values in the Command Window
disp('Gradient Recalled Echo Signal Intensity Values for Patient 1:');

disp('Gray Matter (GM):');
disp(array2table(SI_values_GM, 'VariableNames', {'TR_25', 'TR_50', 'TR_100', 'TR_200'}, 'RowNames', {'Alpha_15', 'Alpha_30', 'Alpha_45', 'Alpha_60', 'Alpha_90'}));

disp('White Matter (WM):');
disp(array2table(SI_values_WM, 'VariableNames', {'TR_25', 'TR_50', 'TR_100', 'TR_200'}, 'RowNames', {'Alpha_15', 'Alpha_30', 'Alpha_45', 'Alpha_60', 'Alpha_90'}));

disp('Cerebrospinal Fluid (CSF):');
disp(array2table(SI_values_CSF, 'VariableNames', {'TR_25', 'TR_50', 'TR_100', 'TR_200'}, 'RowNames', {'Alpha_15', 'Alpha_30', 'Alpha_45', 'Alpha_60', 'Alpha_90'}));
