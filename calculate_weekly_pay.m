%% Calculate weekly Pay Script
disp('This program takes hourly wage rate and hours worked and calculates weekly pay rate');

% Prompt the user to enter hourly wage
promptHourlyWage = 'Enter the hourly wage rate: ';
answerHourlyWage = input(promptHourlyWage);

% Prompt the user to enger hours worked
promptHoursWorked = 'Enter hours worked: ';
answerHoursWorked = input(promptHoursWorked);

weeklyPay = answerHourlyWage * answerHoursWorked;

disp(['The weekly pay is $' num2str(weeklyPay) ]);