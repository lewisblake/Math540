%% Calculate mass consumed script
% This script takes the energy output needed and calculates the mass
% consumed based of E = MC^2

% Define inital quantities: c
speedOfLight_c = 299792458;

% Display program objective to user
disp('This program calculates mass consumed at a nuclear reactor in the course of a year.');

% Get user to input energy output
promptEnergyOutput = 'Enter the energy output (joules): ';
answerEnergyInput = input(promptEnergyOutput);
% Calculate mass consumed
massConsumed = answerEnergyInput/(speedOfLight_c^2);
% Display results
disp(['The ammount of mass consumed is: ', num2str(massConsumed), ' kilograms.']);
