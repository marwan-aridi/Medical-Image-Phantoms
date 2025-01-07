% Assignment 2
% Author: Marwan Aridi, COSC4372
% Note: This Matlab file will generate 4 generated MRIs including legends.
% This file will be useful to generate the phantoms. The file will be
% seperate and is called Q2Phantoms.


% Define compartment parameters
A = [0.9, 0.85, 0.8, 0.7];  % Amount of water in each compartment
T1 = [250, 625, 1000, 1375]; % T1 relaxation times
T2 = [10, 35, 60, 85];       % T2 relaxation times
TE = 10;  % Fixed TE for all images

% Define TR values for each image
TR_values = [50, 250, 1000, 2500];

% Initialize figure for display
figure;

% Loop through TR values and generate images
for img_num = 1:length(TR_values)
    TR = TR_values(img_num);
    
    % Compute the SI for each compartment
    SI = zeros(1, 4);
    for i = 1:4
        SI(i) = A(i) * (1 - exp(-TR / T1(i))) * exp(-TE / T2(i));
    end

    % Create a simple visual representation of the compartments
    compartments = zeros(256, 256);

    % Compartment 1 (center)
    compartments(96:160, 96:160) = SI(1); % Simple square as placeholder

    % Compartment 2 (left)
    compartments(96:160, 32:96) = SI(2);

    % Compartment 3 (right)
    compartments(96:160, 160:224) = SI(3);

    % Compartment 4 (bottom)
    compartments(160:192, 96:160) = SI(4);

    % Display the image in a subplot
    subplot(2, 2, img_num);
    imshow(compartments, []);
    title(sprintf('Image %d: MRI with TR = %d, TE = 10', img_num, TR_values(img_num)));
    legend_text = sprintf('SI1: %.2f, SI2: %.2f, SI3: %.2f, SI4: %.2f', SI(1), SI(2), SI(3), SI(4));
    text(10, 230, legend_text, 'Color', 'white');
end
