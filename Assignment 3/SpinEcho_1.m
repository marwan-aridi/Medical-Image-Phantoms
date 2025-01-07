% Author: Marwan Aridi
% Assignment 3
% This is a Spin Echo file for the first part of the assignment
% Please see README file on how to run this code

% Load patient data and extract slice 90
patient1_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient1.nii');
slice_90_patient1 = patient1_data(:, :, 90);

% Define TR and TE values for Spin Echo table
TR_values = [250, 500, 750, 1000, 2000];   % TR values in ms
TE_values = [20, 40, 60, 80];              % TE values in ms

% Define T1, T2, and P values for each tissue type (GM, WM, CSF)
% Gray Matter (GM)
T1_GM = 1.62; T2_GM = 85.00; P_GM = 105.00;

% White Matter (WM)
T1_WM = 0.92; T2_WM = 80.50; P_WM = 80.00;

% Cerebrospinal Fluid (CSF)
T1_CSF = 10.4; T2_CSF = 1055.00; P_CSF = 150.00;

% Initialize arrays to store signal intensities for each TR and TE combination
SI_values_patient1 = struct('GM', zeros(length(TR_values), length(TE_values)), ...
                            'WM', zeros(length(TR_values), length(TE_values)), ...
                            'CSF', zeros(length(TR_values), length(TE_values)));

% Set up figure for displaying images
figure;
tiledlayout(length(TR_values), length(TE_values), 'TileSpacing', 'compact', 'Padding', 'compact');
title('Patient 1 - Spin Echo Table');

% Calculate signal intensity values for Patient 1 and display images
for i = 1:length(TR_values)
    for j = 1:length(TE_values)
        TR = TR_values(i);
        TE = TE_values(j);
        
        % Initialize the signal intensity map for Patient 1
        SI_map = zeros(size(slice_90_patient1));
        
        % Create tissue masks based on intensity values in the slice
        GM_mask = (slice_90_patient1 == 1);  % Gray Matter
        WM_mask = (slice_90_patient1 == 2);  % White Matter
        CSF_mask = (slice_90_patient1 == 3); % Cerebrospinal Fluid

        % Calculate Signal Intensity for each tissue type
        % Gray Matter
        SI_map(GM_mask) = P_GM * (1 - exp(-TR / (T1_GM * 1000))) * exp(-TE / T2_GM);
        SI_values_patient1.GM(i, j) = mean(SI_map(GM_mask));

        % White Matter
        SI_map(WM_mask) = P_WM * (1 - exp(-TR / (T1_WM * 1000))) * exp(-TE / T2_WM);
        SI_values_patient1.WM(i, j) = mean(SI_map(WM_mask));
        
        % Cerebrospinal Fluid
        SI_map(CSF_mask) = P_CSF * (1 - exp(-TR / (T1_CSF * 1000))) * exp(-TE / T2_CSF);
        SI_values_patient1.CSF(i, j) = mean(SI_map(CSF_mask));

        % Display the image in the corresponding tile
        nexttile((i-1) * length(TE_values) + j);
        imshow(mat2gray(SI_map));  % Normalize and display
        title(sprintf('TR = %d, TE = %d', TR, TE), 'FontSize', 8);
    end
end

% Display results in MATLAB command window
disp('Signal Intensity Values for Patient 1 (Spin Echo):');
disp('Gray Matter (GM):');
disp(SI_values_patient1.GM);
disp('White Matter (WM):');
disp(SI_values_patient1.WM);
disp('Cerebrospinal Fluid (CSF):');
disp(SI_values_patient1.CSF);
