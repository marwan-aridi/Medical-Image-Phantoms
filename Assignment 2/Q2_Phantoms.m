% Assignment 2
% Author: Marwan Aridi, COSC4372
% Note: This Matlab file will generate 4 phantoms that will change based on
% the instructions on Q2.2. 


function generate_mri_images_q2_2_custom_all_in_one()
    % Define the ellipses (same structure as in previous examples)
    ellipses = [
        1, 0.69, 0.92, 0, 0, 0;           % Main large ellipse
       -0.8, 0.6624, 0.874, 0, -0.0184, 0; % Slightly smaller lower ellipse
       -0.2, 0.11, 0.31, 0.3, 0, -18;      % Left smaller ellipse
       -0.2, 0.16, 0.41, -0.3, 0, 18;      % Right smaller ellipse
        0.1, 0.15, 0.18, 0, 0.35, 0;       % Top ellipse
        0.05, 0.046, 0.046, 0.01, -0.1, 0; % Small circle 1
        0.05, 0.046, 0.046, 0.01, -0.3, 0; % Small circle 2
        0.05, 0.023, 0.023, -0.12, -0.58, 0; % Small circle 3
        0.05, 0.023, 0.023, 0.20, -0.58, 0;  % Small circle 4
        0.05, 0.023, 0.023, 0.07, -0.605, 0]; % Small circle 5
    
    % Acquisition parameters from Q2.2
    TR_values = [50, 250, 1000, 2500];  % TR values in ms
    TE_values = [10, 10, 10, 10];       % TE values in ms
    
    % Compartment parameters (A, T1, T2) based on your provided table
    A_values = [0.90, 0.85, 0.80, 0.70];  % A values for the compartments
    T1_values = [250, 625, 1000, 1375];   % T1 values in ms
    T2_values = [10, 35, 60, 85];         % T2 values in ms
    figure;
    
    % Iterate through each set of acquisition parameters and generate MRI images
    for i = 1:4
        TR = TR_values(i);  % Current TR
        TE = TE_values(i);  % Current TE
        
        % Generate the MRI signal intensity map based on the current TR and TE
        [mri_image, SI_values] = generate_signal_intensity_and_si(256, 10, ellipses, A_values, T1_values, T2_values, TR, TE);
        
        % Create a subplot in a 2x2 grid
        subplot(2, 2, i);  % 2 rows, 2 columns, current position
        
        % Display the MRI image
        imagesc(mri_image);    % Use imagesc to scale data and apply colormap
        colormap(gray);        % Use grayscale for MRI image
        colorbar;              % Add a colorbar to show the signal intensity
        title(['MRI Image with TR = ' num2str(TR) ' ms, TE = ' num2str(TE) ' ms']);
        
        % Add a legend to display the acquisition parameters and SI values
        legend_text = sprintf('TR = %d ms, TE = %d ms\nSI values: %.6f, %.6f, %.6f, %.6f', ...
            TR, TE, SI_values(1), SI_values(2), SI_values(3), SI_values(4));
        legend(legend_text, 'Location', 'southoutside');
    end
end

% Function to calculate the MRI signal intensity (SI) based on A, T1, T2, TR, and TE
function [si_map, SI_values] = generate_signal_intensity_and_si(N, M, ellipses, A_values, T1_values, T2_values, TR, TE)
    % Create a blank N x N matrix for the SI-map
    si_map = zeros(N);
    SI_values = zeros(1, 4);  % To store SI for each compartment
    
    % Generate a grid of x-coordinates
    xax = ( (0:N-1)-(N-1)/2 ) / ((N-1)/2); 
    xg = repmat(xax, N, 1);   % X-coordinates
    
    % Iterate through each ellipse and calculate the signal intensity
    for i = 1:M
        A = A_values(mod(i-1, length(A_values)) + 1);   % Amount of water (A)
        T1 = T1_values(mod(i-1, length(T1_values)) + 1); % Longitudinal relaxation time (T1)
        T2 = T2_values(mod(i-1, length(T2_values)) + 1); % Transverse relaxation time (T2)
        
        % Calculate the MRI signal intensity for the current voxel using given TR and TE
        SI = A * (1 - exp(-TR / T1)) * exp(-TE / T2);
        
        % Store the SI value for the current compartment
        if i <= 4
            SI_values(i) = SI;
        end
        
        a = ellipses(i, 2);    % Horizontal semi-axis (a)
        b = ellipses(i, 3);    % Vertical semi-axis (b)
        x0 = ellipses(i, 4);   % X-coordinate of the center
        y0 = ellipses(i, 5);   % Y-coordinate of the center
        phi = ellipses(i, 6);  % Rotation angle (degrees)
        phi = phi * pi / 180;  % Convert degrees to radians
        
        % Shift coordinates to center the ellipse at (x0, y0)
        x = xg - x0;
        y = rot90(xg) - y0;
        
        % Apply the rotation matrix to the coordinates
        cosp = cos(phi); 
        sinp = sin(phi);
        
        % Formula for finding the indices of the pixels inside the ellipse
        idx = find(((x.*cosp + y.*sinp).^2)./a^2 + ((y.*cosp - x.*sinp).^2)./b^2 <= 1);
        
        % Assign the signal intensity to the corresponding pixels
        si_map(idx) = SI;
    end
end

% Call the function to generate MRI images for Q2.2 in a 2x2 layout
generate_mri_images_q2_2_custom_all_in_one();
