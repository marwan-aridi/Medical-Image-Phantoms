% Assignment 2
% Author: Marwan Aridi, COSC4372
% Note: This file will generate 5 output images. The output images are:
% A-map
% T1-map (greyscale image)
% T2-map (greyscale image)
% T1-map (color)
% T2-map (color)




% This is a combined script for A-map, T1-map, and T2-map (grayscale and colored)

% Main script to generate and display all maps for the modified Shepp-Logan phantom
function generate_all_maps()
    % Define the parameters for the ellipses, matching each to the compartment
    ellipses = [
        1, 0.69, 0.92, 0, 0, 0;           % Main ellipse
       -0.8, 0.6624, 0.874, 0, -0.0184, 0; % Compartment 1 
       -0.2, 0.11, 0.31, 0.3, 0, -18;      % Compartment 2 
       -0.2, 0.16, 0.41, -0.3, 0, 18;      % Compartment 2 
        0.1, 0.15, 0.18, 0, 0.35, 0;       % Compartment 3 
        0.05, 0.046, 0.046, 0.01, -0.1, 0; % Small circle 3
        0.05, 0.046, 0.046, 0.01, -0.3, 0; % Small circle 3
        0.05, 0.023, 0.023, -0.12, -0.58, 0; % Small circle 4
        0.05, 0.023, 0.023, 0.20, -0.58, 0;  % Small circle 4
        0.05, 0.023, 0.023, 0.07, -0.605, 0]; % Small circle 4

    compartments = [1, 2, 3, 3, 4, 4, 4, 4, 4, 4];  
    T1_values = 250 + (compartments - 1) * 375; %Formula for T1
    T2_values = 10 + (compartments - 1) * 25;   %Formula for T2
    A_values = [0.9, 0.8, 0.6, 0.6, 0.9, 0.9, 0.9, 0.4, 0.4, 0.4];  % A-map values

    % Generate the A-map
    A_map = phantom1_A(256, 10, ellipses, A_values);
    figure, imshow(A_map, []);  % Display the A-map
    title('A-map');  % Add the title "A-map"

    % Generate the T1-map with colors
    T1_map_colored = phantom1_T1(256, 10, ellipses, T1_values);
    figure, imagesc(T1_map_colored);  % Use imagesc to display with colors
    colormap(jet);                    % Apply the 'jet' colormap for vibrant colors
    colorbar;                         % Display color bar for reference
    title('T1-map with Colors');

    % Generate the T2-map with colors
    T2_map_colored = phantom1_T2(256, 10, ellipses, T2_values);
    figure, imagesc(T2_map_colored);  % Use imagesc to display with colors
    colormap(hot);                    % Apply the 'hot' colormap for distinct colors
    colorbar;                         % Display color bar for reference
    title('T2-map with Colors');

    % Generate the grayscale T1-map
    T1_map = phantom1_T1(256, 10, ellipses, T1_values);
    figure, imshow(T1_map, []);  % Display the T1-map
    title('T1-map');  % Add the title "T1-map"

    % Generate the grayscale T2-map
    T2_map = phantom1_T2(256, 10, ellipses, T2_values);
    figure, imshow(T2_map, []);  % Display the T2-map
    title('T2-map');  % Add the title "T2-map"
end

% Modified phantom1 function to create the A-map
function A_map = phantom1_A(N, M, ellipses, A_values)
    A_map = zeros(N);
    xax = ( (0:N-1)-(N-1)/2 ) / ((N-1)/2); 
    xg = repmat(xax, N, 1);
    
    for i = 1:M
        A = A_values(i);  % Assign value of A for the compartment
        a = ellipses(i, 2);
        b = ellipses(i, 3);
        x0 = ellipses(i, 4);
        y0 = ellipses(i, 5);
        phi = ellipses(i, 6);
        phi = phi * pi / 180;
        
        x = xg - x0;
        y = rot90(xg) - y0;
        
        cosp = cos(phi); 
        sinp = sin(phi);
        
        idx = find(((x.*cosp + y.*sinp).^2)./a^2 + ((y.*cosp - x.*sinp).^2)./b^2 <= 1);
        A_map(idx) = A;
    end
end

% Modified phantom1 function to create the T1-map
function T1_map = phantom1_T1(N, M, ellipses, T1_values)
    T1_map = zeros(N);
    xax = ( (0:N-1)-(N-1)/2 ) / ((N-1)/2); 
    xg = repmat(xax, N, 1);
    
    for i = 1:M
        T1 = T1_values(i);
        a = ellipses(i, 2);
        b = ellipses(i, 3);
        x0 = ellipses(i, 4);
        y0 = ellipses(i, 5);
        phi = ellipses(i, 6);
        phi = phi * pi / 180;
        
        x = xg - x0;
        y = rot90(xg) - y0;
        
        cosp = cos(phi); 
        sinp = sin(phi);
        
        idx = find(((x.*cosp + y.*sinp).^2)./a^2 + ((y.*cosp - x.*sinp).^2)./b^2 <= 1);
        T1_map(idx) = T1;
    end
end

% Modified phantom1 function to create the T2-map
function T2_map = phantom1_T2(N, M, ellipses, T2_values)
    T2_map = zeros(N);
    xax = ( (0:N-1)-(N-1)/2 ) / ((N-1)/2); 
    xg = repmat(xax, N, 1);
    
    for i = 1:M
        T2 = T2_values(i);
        a = ellipses(i, 2);
        b = ellipses(i, 3);
        x0 = ellipses(i, 4);
        y0 = ellipses(i, 5);
        phi = ellipses(i, 6);
        phi = phi * pi / 180;
        
        x = xg - x0;
        y = rot90(xg) - y0;
        
        cosp = cos(phi); 
        sinp = sin(phi);
        
        idx = find(((x.*cosp + y.*sinp).^2)./a^2 + ((y.*cosp - x.*sinp).^2)./b^2 <= 1);
        T2_map(idx) = T2;
    end
end

% Call the function to generate and display all maps
generate_all_maps();
