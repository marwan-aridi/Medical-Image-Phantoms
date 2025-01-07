% Assignment 2
% Author: Marwan Aridi, COSC4372
% Note: This Matlab file will generate 2 plots for questions Q2.3 and Q2.4



% This is the main script to generate two separate output graphs
function generate_two_graphs()
    % First Graph: Signal Intensity (SI) vs Echo Time (TE) for Compartment 2
    generate_si_vs_te();
    
    % Second Graph: Signal Intensity (SI) vs Repetition Time (TR) for all compartments
    generate_si_vs_tr();
end



% Function to generate the first graph: SI vs TE for Compartment 2
function generate_si_vs_te()
    % Given values for Compartment 2
    A_compartment_2 = 0.85;  % A value for Compartment 2
    T1_compartment_2 = 625;  % T1 value for Compartment 2 in ms
    T2_compartment_2 = 35;   % T2 value for Compartment 2 in ms

    % TR fixed at 250 ms
    TR_fixed = 250;

    % TE values for the question
    TE_values = [10, 40, 80, 150];  % TE values in ms

    % Preallocate array to store SI values for each TE
    SI_values = zeros(1, length(TE_values));

    % Calculate the SI for each TE value
    for i = 1:length(TE_values)
        TE = TE_values(i);
        SI_values(i) = A_compartment_2 * (1 - exp(-TR_fixed / T1_compartment_2)) * exp(-TE / T2_compartment_2);
    end

    % Plotting SI vs TE
    figure;
    plot(TE_values, SI_values, '-o', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('Echo Time (TE) [ms]');
    ylabel('Signal Intensity (SI)');
    title('Signal Intensity vs Echo Time (TE) for Compartment 2 (TR = 250 ms)');
    grid on;
end

% Function to generate the second graph: SI vs TR for all compartments
function generate_si_vs_tr()
    % SI values from the table for each compartment and TR value
    SI = [
        0.0600167, 0.2092897, 0.3250273, 0.3310765;  % Compartment 1
        0.0491099, 0.2105849, 0.5097931, 0.6270565;  % Compartment 2
        0.0330267, 0.1497929, 0.4280628, 0.6215986;  % Compartment 3
        0.0222228, 0.1034567, 0.3215926, 0.5212936   % Compartment 4
    ];

    % TR values
    TR_values = [50, 250, 1000, 2500];

    % Plotting the SI vs TR for each compartment
    figure;
    hold on;
    plot(TR_values, SI(1, :), '-o', 'DisplayName', 'Compartment 1');
    plot(TR_values, SI(2, :), '-s', 'DisplayName', 'Compartment 2');
    plot(TR_values, SI(3, :), '-d', 'DisplayName', 'Compartment 3');
    plot(TR_values, SI(4, :), '-^', 'DisplayName', 'Compartment 4');
    hold off;

    % Labels and title
    xlabel('Repetition Time (TR) [ms]');
    ylabel('Signal Intensity (SI)');
    title('Signal Intensity vs Repetition Time for Different Compartments');
    legend('show');
    grid on;
end

% Call the main function to generate the graphs
generate_two_graphs();
