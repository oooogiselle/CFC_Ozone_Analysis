%% Preliminary Questions
% What is the total amount of ozone and CHC stored in the system at the beginning of the simulation?

% The total amount of ozone at the beginning is 3136e + 09 kg. The total amount of CHC stored in the system is 2kg at the beginning. 

% How would they evolve over time in the absence of CFC emissions (Em = 0
% GtC/yr)? Describe briefly how you expect them to evolve over time in this case.

% Without CFC emissions, the CFC in troposphere and the stratosphere will transfer between each other. The photochemical loss will gradually lead to a decay in CFCs in both spheres eventually. 
% Without CFC emissions, the reduction of CFC (as discussed above) will lead to a decreased rate in ozone depletion, which would potentially allow the ozone levels to recover and stablize over time. 

% How would they evolve over time in the presence of large CFC emissions? Describe
% briefly how you expect them to evolve over time in this case

% In the presence of large CFC emissions, the CFCs in the troposphere will approximately increase linearly with a rate of Em that will eventually stabilize. 
% Initially, CF Cstrat will also start increasing due to the transfer from the troposphere. 
% As CF Cstrat increases, the rate of ozone loss will increase, leading to a more rapid decrease in ozone levels which will eventually stabilize.

%% Part 1
%declare constants
kst = 0.25; %yr-1
kts = 0.05; %yr-1
ko = 3.259; %yr-1
kco = 8.33 * 10^-10; %yr-1
t_photo = 20; %yr
em = 111 * 10^6; %kg yr-1
Op = 28 * 10^9 * 365; %yr-1

% initial condition 
CFC_trop_initial = 1; %kg
CFC_strat_initial = 1; %kg
O3_intiial = 3136 *10^9; %kg

%time parameters
dt = 0.01; %year
t_final = 1000; %year
n_steps = t_final/dt;
time = linspace(0, t_final, n_steps); %declare an array for time 

%initialize arrays with initial conditions
Ctrop_array = CFC_trop_initial;
Cstrat_array = CFC_strat_initial;
O3_array = O3_intiial;

%Euler's Method
for k = 2:n_steps
    %system of eqs
    dCFC_trop = em + kst * CFC_strat_initial - kts * CFC_trop_initial;
    dCFC_strat = kts * CFC_trop_initial - kst * CFC_strat_initial - CFC_strat_initial / t_photo;
    dO3 = Op - ko * O3_intiial - kco * O3_intiial * CFC_strat_initial;

    %step-wise increase
    CFC_trop_initial = CFC_trop_initial + dCFC_trop * dt;
    CFC_strat_initial = CFC_strat_initial + dCFC_strat * dt;
    O3_intiial = O3_intiial + dO3 * dt;

    %assigning value to the array
    Ctrop_array(k) = CFC_trop_initial;
    Cstrat_array(k) = CFC_strat_initial;
    O3_array(k) = O3_intiial;
end

%plots using subplots!
figure;
subplot(1,3,1);
plot(time, Ctrop_array, 'LineWidth', 2);
title('time vs CFC trop');
xlabel('time (years)');
ylabel('CFC in trop (kg)');
legend ('CFC trop');
grid on;

subplot(1,3,2);
plot(time, Cstrat_array, 'LineWidth', 2);
title('time vs CFC strat');
xlabel('time (years)');
ylabel('CFC in strat (kg)');
legend ('CFC strat');
grid on;

subplot(1,3,3);
plot(time, O3_array, 'LineWidth', 2);
title('time vs ozone');
xlabel('time (years)');
ylabel('ozone (kg)');
legend ('ozone');
grid on;

% extra credit
figure;
yyaxis left;
plot(time, Ctrop_array,'g-', time, Cstrat_array, 'r-');
ylabel('CFC level (kg)');

yyaxis right;
plot(time, O3_array, 'b-');
ylabel('Ozone (kg)');

xlabel('time (years)');
title('time vs CFCs in trop, strat, and Ozone (extra credit)');
legend('CFC in trop', 'CFC in strat','Ozone');
grid on;

% It takes about 600 years to reach a steady state (from the plot). 
% The initial value of CFC in the troposphere is 1 kg. The steady-state value appears to be around 1.3e+10kg.
% The initial value of CFC in the stratosphere is 1kg. The steady-state value is about 2.25e+9 kg. 
% The initial value of ozone is 3.136e+12 kg. The steady-state value seems to be around 2e+12.
% The final amount of CFC in the troposphere is about 1.3e+10kg and the final amount of steady-state value of CFC in the stratosphere is 2.25e+9 kg.
% The total emission over this period is 1.11e+11 kg (111e+6 kg/yr * e+3 yr).
% The final amount of CFC is greater than the CFC emissions over that period. 
%% Part 2 

% initial condition 
CFC_trop_initial = 1; 
CFC_strat_initial = 1;
O3_initial = 3136 *10^9;
em = 111 * 10^6;

% call function that calculates the first 100 years with emissions 
[Ctrop_array_withE,Cstrat_array_withE, O3_array_withE, time1] = Emissions100(CFC_trop_initial, CFC_strat_initial, O3_initial, 0, 100, em); %call the function to calculate 100 years with emissions
[Ctrop_array_withoutE,Cstrat_array_withoutE, O3_array_withoutE, time2] = Emissions100(Ctrop_array_withE(end), Cstrat_array_withE(end), O3_array_withE(end), 100, 200, 0);

% Combine time arrays and values
combined_time = [time1, time2(2:end)]; %start from 2 to prevent overlap 
combined_Ctrop = [Ctrop_array_withE, Ctrop_array_withoutE(2:end)];
combined_Cstrat = [Cstrat_array_withE, Cstrat_array_withoutE(2:end)];
combined_O3 = [O3_array_withE, O3_array_withoutE(2:end)];

% Plot the results
figure;
subplot(3,1,1);
plot(combined_time, combined_Ctrop);
title('CFC in trop before and after the ban');
xlabel('time (years)');
ylabel('CFC in trop (kg)');
grid on;

subplot(3,1,2);
plot(combined_time,combined_Cstrat);
title('CFC in strat before and after the ban');
xlabel('time (years)');
ylabel('CFC in strat (kg)');
grid on;

subplot(3,1,3);
plot(combined_time, combined_O3);
title('ozone before and after the ban');
xlabel('time (years)');
ylabel('ozone (kg)');
grid on;

% extra credit
figure;
yyaxis left;
plot(combined_time, combined_Ctrop, combined_time, combined_Cstrat);
ylabel('CFC level (kg)');

yyaxis right;
plot(combined_time, combined_O3);
ylabel('Ozone level (kg)');

title('CFC in trop, strat, and ozone before and after ban (extra credit)');
xlabel('time (years)');
legend('CFC in trop', 'CFC in strat','Ozone');
grid on;

% Ozone did not recover to its pre-CFCs value, which is 3.136e+12 kg.
% The final ozone level after 100 years of ban is 2.74e+12kg.
% There is still 3.96e+11kg of ozone depletion remaining after 100 years of no CFC emission. 

%% part 3
%load data
ozone_data = load('Ozone.txt');

% plot data points loaded
figure;
plot(ozone_data(:,1), ozone_data(:,2), 'o');
hold on;

% initial condition 
CFC_trop_initial = 1; 
CFC_strat_initial = 1;
O3_initial = 3136 *10^9;
em = 111 * 10^6;

% call functions to get the level 
[Ctrop_array_withE,Cstrat_array_withE, O3_array_withE, ~] = Emissions100(CFC_trop_initial, CFC_strat_initial, O3_initial, 0, 100, em); %call the function to calculate 100 years with emissions
[~,~, O3_array_withoutE, time2] = Emissions100(Ctrop_array_withE(end), Cstrat_array_withE(end), O3_array_withE(end), 100, 200, 0);

% Extract data between years 20 and 30
start_index = find(time2 >= 120, 1); % Find the index corresponding to year 120
end_index = find(time2 <= 130, 1, 'last'); % Find the index corresponding to year 130

time_range = time2(start_index:end_index) - 100; % Adjust to get the correct time range (20 to 30 years)
O3_range = O3_array_withoutE(start_index:end_index);

% Plot the simulated ozone levels for the specified range
plot(time_range, O3_range);
xlabel('Time (years)');
ylabel('O3 (kg)');
title('Ozone Levels Between 120 and 130 Years');

% run bisection method to find the least squared fit 
a = 0;
b = em; % em max 
tol = 0.01;

while (abs(b-a)>0.01)
    c = (a+b)/2;

    error_a = leastSquares(a,ozone_data, Ctrop_array_withE, Cstrat_array_withE, O3_array_withE);
    error_b = leastSquares(b,ozone_data, Ctrop_array_withE, Cstrat_array_withE, O3_array_withE);
    error_c = leastSquares(c,ozone_data, Ctrop_array_withE, Cstrat_array_withE, O3_array_withE);

    if error_a >= error_b
        a = c;
    else
        b = c;
    end
end

optimized_em = c; %optimized emission 
least_error = error_c; % least squared error 
 
% call functions to get the level with the optimized emission 
[~,~, O3_leastsquares, time4] = Emissions100(Ctrop_array_withE(end), Cstrat_array_withE(end), O3_array_withE(end), 100, 200, optimized_em);

% Extract data between years 20 and 30
start_idx = find(time4 >= 120, 1); % Find the index corresponding to year 120
end_idx = find(time4 <= 130, 1, 'last'); % Find the index corresponding to year 130

time_range_ls = time4(start_idx:end_idx) - 100; % Adjust to get the correct time range (20 to 30 years)
O3_range_ls = O3_leastsquares(start_idx:end_idx);

% Plot the least squares fit
plot(time_range_ls, O3_range_ls);
legend('Measured Data', 'Simulated Data', 'Least Squares Fit');
hold off;

fprintf('The emission that gives you the least squares fit is %.2d kg/year with an error of %.2d \n',optimized_em, least_error);

%% Extra Credit (keeping track of chemical loss and ozone loss)

% initial condition 
CFC_trop_initial = 1; 
CFC_strat_initial = 1;
O3_initial = 3136 *10^9;
em = 111 * 10^6;
t_photo = 20; %yr

% call function that calculates the first 100 years with emissions 
[Ctrop_array_withE,Cstrat_array_withE, O3_array_withE, time1] = Emissions100(CFC_trop_initial, CFC_strat_initial, O3_initial, 0, 100, em); %call the function to calculate 100 years with emissions
[Ctrop_array_withoutE,Cstrat_array_withoutE, O3_array_withoutE, time2] = Emissions100(Ctrop_array_withE(end), Cstrat_array_withE(end), O3_array_withE(end), 100, 200, 0);

% Combine time arrays and values
combined_time = [time1, time2(2:end)]; %start from 2 to prevent overlap 

% Calculate chemical loss
chemical_loss_withE = Cstrat_array_withE ./ t_photo;
chemical_loss_withoutE = Cstrat_array_withoutE ./ t_photo;
combined_chemical_loss = [chemical_loss_withE, chemical_loss_withoutE(2:end)];

% Calculate ozone loss from CFC chemistry
kco = 8.33 * 10^-10; % Rate constant for ozone destruction
ozone_loss_withE = kco * O3_array_withE .* Cstrat_array_withE;
ozone_loss_withoutE = kco * O3_array_withoutE .* Cstrat_array_withoutE;
combined_ozone_loss = [ozone_loss_withE, ozone_loss_withoutE(2:end)];

% Plot chemical loss
figure;
subplot(2, 1, 1);
plot(combined_time, combined_chemical_loss, 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Chemical Loss of CFCs');
title('Chemical Loss of CFCs in the Stratosphere Over Time');
legend('Chemical Loss of CFCs');

% Plot ozone loss
subplot(2, 1, 2);
plot(combined_time, combined_ozone_loss, 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Ozone Loss');
title('Ozone Loss from CFC Chemistry Over Time');
legend('Ozone Loss');

% Describe the results
disp('During the first 100 years with emissions, both the chemical loss of CFCs and ozone loss increase.');
disp('After the ban on CFC emissions, both the chemical loss of CFCs and ozone loss gradually decrease as the concentration of CFCs in the stratosphere declines.');

%% ODE solver
% [t,y] = ode45(odefun,tspan,y0), where tspan = [t0 tf], integrates the system of differential equations y′=f(t,y) from t0 to tf with initial conditions y0. 
% Each row in the solution array y corresponds to a value returned in column vector t.

%declare constants
kst = 0.25; %yr-1
kts = 0.05; %yr-1
ko = 3.259; %yr-1
kco = 8.33 * 10^-10; %yr-1
t_photo = 20; %yr
em = 111 * 10^6; %yr-1
Op = 28 * 10^9 * 365; %yr-1

%initial conditions
y0 = [1,1,3136 *10^9];

% declare time array
t0 = 0;
tf = 1000;
tspan = [t0 tf];

% define system of functions
odefun = @(t,y) [
    em + kst * y(2) - kts * y(1);                    % dCFC_trop/dt
    kts * y(1) - kst * y(2) - y(2) / t_photo;        % dCFC_strat/dt
    Op - ko * y(3) - kco * y(3) * y(2)               % dO3/dt
];

% call built-in ODE function
[t,y] = ode45(odefun, tspan,y0);

% Extract solutions 
Ctrop_solution = y(:,1);
Cstrat_solution = y(:,2);
O3_solution = y(:,3);

% Plot! 
figure;
yyaxis left;
plot(t,Ctrop_solution,t,Cstrat_solution);
ylabel('CFC level (kg)');

yyaxis right;
plot(t,O3_solution);
ylabel('Ozone (kg)');

xlabel('time (years)');
title('time vs CFC and Ozone');



