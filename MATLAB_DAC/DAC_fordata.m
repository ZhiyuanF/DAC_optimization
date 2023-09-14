function [results] = DAC_fordata(price_cap, price_data, X, k, parameters)

    % Fixed parameters
    X_hat = parameters(1,1); % Paying S dollar per each cycle for material/operation cost
    S = parameters(1,2); % Paying S dollar per each cycle for material/operation cost
    pi_co2 = parameters(1,3); % Selling co2 at pi_co2 per ton, incentive
    P_a = parameters(1,6); % Correction of absorption power consumption by the installed capacity
    P_d = parameters(1,7); % Correction of desorption power consumption by the installed capacity
    beta_a_1 = parameters(1,8); % Piecewise linear approximation for absorption, first-order term
    beta_a_2 = parameters(1,9); % Piecewise linear approximation for absorption, second-order term
    beta_d_1 = parameters(1,10); % Piecewise linear approximation for desorption, first-order term
    beta_d_2 = parameters(1,11); % Piecewise linear approximation for desorption, second-order term

    % Initialize variables
    results = cell(numel(price_data),9);

    for t = 1:numel(price_data)
        if price_data(t) > price_cap %if price is higher than lambda, idle
            u = false;
            v = false;
            z = false;

        elseif X > 0.8*X_hat %if saturation is high, desorb
            u = false;
            v = true;
            z = false;
            k = false;

        elseif X < 0.05*X_hat %if saturation is low, adsorb

            % check if most recently desorbing,
            % and if so, indicate a new cycle has begun
            u = true;
            v = false;
            z = ~k;
            k = true;
        else
            % continue adsorbing or desorbing as previously
            u = k;
            v = ~k;
            z = false;
        end

        % calculating adsorption and desorption rates, state of saturation,
        % tons of CO2 sold, and profit
        a = (beta_a_1 + beta_a_2 * X) * u;
        d = (beta_d_1 + beta_d_2 * X) * v;
        X = X + a - d;
        profit = (pi_co2 * d) - (price_data(t) * (P_a * u + P_d * v) + S * z);

        % Store the values in the results array
        results(t, :) = {u, v, k, z, price_data(t), X, a, d, profit};
    end
end