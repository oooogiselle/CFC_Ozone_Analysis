function [Ctrop_array,Cstrat_array, O3_array, time] = Emissions100(CFC_trop_initial,CFC_strat_initial, O3_initial, t_initial, t_final, em)
%EMISSIONS100 calculates the CFC in troposphere and stratosphere
%and ozone level after certain number of  years of emission. 
%   This function simulates the impact of continuous emissions on the
%   concentrations of CFCs in troposphere and stratosphere, as well as the resulting changes in ozone levels over 
%   a specified time period using Euler's method for numerical integration.

%declare constants
kst = 0.25;
kts = 0.05;
ko = 3.259;
kco = 8.33 * 10^-10;
t_photo = 20;
Op = 28 * 10^9 * 365;


%time parameters
dt = 0.01; %year
n_steps = (t_final-t_initial)/dt;
time = linspace(t_initial, t_final, n_steps); %declare an array for time 

% Initialize arrays with initial conditions
Ctrop_array = zeros(1, n_steps);
Cstrat_array = zeros(1, n_steps);
O3_array = zeros(1, n_steps);

Ctrop_array(1) = CFC_trop_initial;
Cstrat_array(1) = CFC_strat_initial;
O3_array(1) = O3_initial;

%Euler's Method
for k = 2:n_steps
    % System of equations
    dCFC_trop = em + kst * Cstrat_array(k-1) - kts * Ctrop_array(k-1);
    dCFC_strat = kts * Ctrop_array(k-1) - kst * Cstrat_array(k-1) - Cstrat_array(k-1) / t_photo;
    dO3 = Op - ko * O3_array(k-1) - kco * O3_array(k-1) * Cstrat_array(k-1);
    
    % Step-wise increase
    Ctrop_array(k) = Ctrop_array(k-1) + dCFC_trop * dt;
    Cstrat_array(k) = Cstrat_array(k-1) + dCFC_strat * dt;
    O3_array(k) = O3_array(k-1) + dO3 * dt;
end
end


