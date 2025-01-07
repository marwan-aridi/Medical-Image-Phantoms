% Author: Marwan Aridi
% Assignment 3
% This is a T1 Inversion file for the first part of the assignment
% Please see README file on how to run this code


% Load patient data and extract slice 90
patient1_data = niftiread('/Users/marwanaridi/Desktop/COSC4372-6370-Assignment_Synthesis_ML-OASIS/Data/patient1.nii');
slice_90_patient1 = patient1_data(:, :, 90);

% Define TR and TI values for T1 Inversion Recovery table
TR_values = [1000, 2000];  % TR values in ms
TI_values = [50, 100, 250, 500, 750];  % TI values in ms

% Define T1, T2, and P values for each tissue type (GM, WM, CSF)
% Gray Matter (GM)
T1_GM = 1.62; T2_GM = 85.00; P_GM = 105.00;

% White Matter (WM)
T1_WM = 0.92; T2_WM = 80.50; P_WM = 80.00;

% Cerebrospinal Fluid (CSF)
T1_CSF = 10.4; T2_CSF = 1055.00; P_CSF = 150.00;

% Convert T1 values to milliseconds for consistency with TR and TI
T1_GM = T1_GM * 1000; % Convert GM T1 to ms
T1_WM = T1_WM * 1000; % Convert WM T1 to ms
T1_CSF = T1_CSF * 1000; % Convert CSF T1 to ms

% Initialize matrices to store results for Patient 1
gm_values_patient1 = zeros(length(TR_values), length(TI_values));
wm_values_patient1 = zeros(length(TR_values), length(TI_values));
csf_values_patient1 = zeros(length(TR_values), length(TI_values));

% Set up figure for displaying images
figure;
tiledlayout(length(TI_values), length(TR_values), 'TileSpacing', 'compact', 'Padding', 'compact');
title('Patient 1 - T1 Inversion Recovery Table');

% Calculate signal intensities for each combination of TR and TI for Patient 1
for i = 1:length(TI_values)
    for j = 1:length(TR_values)
        TR = TR_values(j);
        TI = TI_values(i);
        
        % Initialize the signal intensity map for Patient 1
        SI_map = zeros(size(slice_90_patient1));
        
        % Create tissue masks based on intensity values in the slice
        GM_mask = (slice_90_patient1 == 1);  % Gray Matter
        WM_mask = (slice_90_patient1 == 2);  % White Matter
        CSF_mask = (slice_90_patient1 == 3); % Cerebrospinal Fluid

        % Calculate Signal Intensity for each tissue type
        % Gray Matter
        SI_map(GM_mask) = P_GM * (1 - 2 * exp(-TI / T1_GM) + exp(-TR / T1_GM));
        gm_values_patient1(j, i) = mean(SI_map(GM_mask));
        
        % White Matter
        SI_map(WM_mask) = P_WM * (1 - 2 * exp(-TI / T1_WM) + exp(-TR / T1_WM));
        wm_values_patient1(j, i) = mean(SI_map(WM_mask));
        
        % Cerebrospinal Fluid
        SI_map(CSF_mask) = P_CSF * (1 - 2 * exp(-TI / T1_CSF) + exp(-TR / T1_CSF));
        csf_values_patient1(j, i) = mean(SI_map(CSF_mask));
        
        % Display the image in the corresponding tile
        nexttile((i-1) * length(TR_values) + j);
        imshow(mat2gray(SI_map));  % Normalize and display
        title(sprintf('TR = %d, TI = %d', TR, TI), 'FontSize', 8);
    end
end

% Display results in MATLAB command window
disp('Gray Matter (GM) Signal Intensities for Patient 1:');
disp(gm_values_patient1);

disp('White Matter (WM) Signal Intensities for Patient 1:');
disp(wm_values_patient1);

disp('Cerebrospinal Fluid (CSF) Signal Intensities for Patient 1:');
disp(csf_values_patient1);
