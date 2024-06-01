function error = leastSquares(em, ozone_data, Ctrop_array_withE, Cstrat_array_withE, O3_array_withE)
%LEASTSQUARES calculates the least square
%   em: emission value
%   ozone_data: measured ozone data
%   Ctrop_array_withE, Cstrat_array_withE, O3_array_withE: arrays from previous emission calculations

% run the Euler's method for the input em
[~,~, O3sim, ~] = Emissions100(Ctrop_array_withE(end), Cstrat_array_withE(end), O3_array_withE(end), 100, 200, em);

% allocate an error array
errorarray = zeros(1,length(ozone_data));

% for loop that populates the error array
for i = 1:length(ozone_data)
    temperror = ozone_data(i,2) - (O3sim(ozone_data(i,1)*100));
    errorarray(i) = temperror;
end

% find the sum of least squares
error = sum((errorarray).^2);

end

