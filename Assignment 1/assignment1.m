% Assignment 1
% Author: Marwan Aridi, COSC4372
% Note: This is a one MATLAB file that will generate all 4 images that are
% required for this assignment. Please see README file for details on how
% to run it.


% Main script: Run all functions to generate the four output pictures
image1(); % This image is for Q1.1
image2(); % This image is for Q2.1
image3(); % This image is for Q2.2
image4(); % This image is for Q2.3


% Main script to generate and display all four images

% Image 1 function
function image1()
    % Define the parameters for the ellipses
    % Each row specifies: [Intensity, a, b, x0, y0, phi]
    ellipses = [1, 0.69, 0.92, 0, 0, 0;   
               -0.8, 0.6624, 0.8740, 0, -0.0184, 0;  
               -0.2, 0.11, 0.31, 0.22, 0, -18;       
               -0.2, 0.16, 0.41, -0.22, 0, 18;       
                0.1, 0.21, 0.25, 0, 0.35, 0;         
                0.1, 0.046, 0.046, 0, 0.1, 0;        
                0.1, 0.046, 0.046, 0, -0.1, 0;       
                0.1, 0.046, 0.023, -0.08, -0.605, 0; 
                0.1, 0.023, 0.023, 0, -0.606, 0;    
                0.1, 0.023, 0.046, 0.06, -0.605, 0]; 

    % Call phantom1 to create and display the image with 10 ellipses
    custom_phantom = phantom1(256, 10, ellipses);
    figure, imshow(custom_phantom, []);
end

% Image 2 function
function image2()
    % Define the parameters for the ellipses with edited overlap
    ellipses = [
        1, 0.69, 0.92, 0, 0, 0;           
       -0.8, 0.6624, 0.874, 0, -0.0184, 0; 
       -0.2, 0.11, 0.31, 0.3, 0, -18;      
       -0.2, 0.16, 0.41, -0.3, 0, 18;      
        0.1, 0.15, 0.18, 0, 0.35, 0;       
        0.05, 0.046, 0.046, 0.01, -0.1, 0; 
        0.05, 0.046, 0.046, 0.01, -0.3, 0; 
        0.05, 0.023, 0.023, -0.12, -0.58, 0; 
        0.05, 0.023, 0.023, 0.20, -0.58, 0;  
        0.05, 0.023, 0.023, 0.07, -0.605, 0]; 
    
    % Call phantom1 to create and display the image without overlap
    custom_phantom = phantom1(256, 10, ellipses);
    figure, imshow(custom_phantom, []);
end

% Image 3 function
function image3()
    % Define parameters including two circles outside the main structure
    ellipses = [
        1, 0.69, 0.92, 0, 0, 0;           
       -0.8, 0.6624, 0.874, 0, -0.0184, 0; 
       -0.2, 0.11, 0.31, 0.3, 0, -18;      
       -0.2, 0.16, 0.41, -0.3, 0, 18;      
        0.1, 0.15, 0.18, 0, 0.35, 0;       
        0.05, 0.046, 0.046, 0, -0.1, 0;    
        0.05, 0.046, 0.046, 0, -0.3, 0;    
        0.05, 0.023, 0.023, -0.12, -0.58, 0; 
        0.05, 0.023, 0.023, 0.13, -0.58, 0;  
        0.05, 0.023, 0.023, 0.06, -0.605, 0; 
        0.9, 0.1, 0.1, -0.88, 0.78, 0;   % Circle 1    
        0.9, 0.1, 0.1, 0.88, 0.78, 0];   % Circle 2   
    
    % Call phantom1 to create and display the image with 12 ellipses
    custom_phantom = phantom1(256, 12, ellipses);
    figure, imshow(custom_phantom, []);
end

% Image 4 function
function image4()
    % Define parameters, removing Q2.2 shapes and adding three concentric circles
    ellipses = [
        1, 0.69, 0.92, 0, 0, 0;           % Main large ellipse
       -0.8, 0.6624, 0.8740, 0, -0.0184, 0; % Slightly smaller lower ellipse
        0.357, 0.1, 0.1, 0, 0, 0;         % Inner concentric circle 1
        0.257, 0.2, 0.2, 0, 0, 0;         % Inner concentric circle 2
        0.157, 0.3, 0.3, 0, 0, 0];        % Outer concentric circle
    
    % Call phantom1 to create and display the image with 5 ellipses
    custom_phantom = phantom1(256, 5, ellipses);
    figure, imshow(custom_phantom, []);
end



% The source code below was taken from MATLAB with added edits
% Function to create the phantom used by all 4 images
function custom_phantom = phantom1(N, M, ellipses)
    % Create a blank N x N matrix for the phantom image
    custom_phantom = zeros(N);
    
    % Generate a grid of x-coordinates
    xax = ( (0:N-1)-(N-1)/2 ) / ((N-1)/2); 
    xg = repmat(xax, N, 1);   % X-coordinates
    
    % Iterate through each ellipse and incorporate it into the phantom
    for i = 1:M
        A = ellipses(i, 1);    % Intensity of the ellipse
        a = ellipses(i, 2);    % Horizontal semi-axis (a)
        b = ellipses(i, 3);    % Vertical semi-axis (b)
        x0 = ellipses(i, 4);   % X-coordinate of the center
        y0 = ellipses(i, 5);   % Y-coordinate of the center
        phi = ellipses(i, 6);  % Rotation angle (degrees)
        phi = phi * pi / 180;  % Convert degrees to radians
        
        % Shift coordinates to center the ellipse at (x0, y0)
        x = xg - x0;
        y = rot90(xg) - y0;
        
        % Applying the rotation matrix to the coordinates
        cosp = cos(phi); 
        sinp = sin(phi);
        
        % formula for finding the indices of the pixels inside the ellipse
        idx = find(((x.*cosp + y.*sinp).^2)./a^2 + ((y.*cosp - x.*sinp).^2)./b^2 <= 1);
        
        % formula for adding the ellipse's intensity to the corresponding pixels
        custom_phantom(idx) = custom_phantom(idx) + A;
    end
end
