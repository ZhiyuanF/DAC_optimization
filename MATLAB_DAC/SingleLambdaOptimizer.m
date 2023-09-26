function [best_lambda_opt, best_prof] = SingleLambdaOptimizer(price_data, X, k, parameters)

    % define the function that wraps the DAC_Cycle function
    objectiveFunc = @(lambda)(-DAC_foropt(lambda, price_data, X, k, parameters));
    
    % initialize variables
    best_prof = 0; % initialize with zero
    best_lambda_opt = 0; % initialize with zero
    counter = 0; % initialize with zero
    
    for lambda_guess = -10:1:60
                
        % perform optimization using fminsearch
        [lambda_opt, fake_profit] = fminsearch(objectiveFunc, lambda_guess);
    
        % update the best optimized output and corresponding lambda_opt
        if -fake_profit > best_prof
            best_prof = -fake_profit;
            best_lambda_opt = lambda_opt;
            counter = counter + 1;
        end
    end

    %double check if the best lambda opt has ever been changed, if not, use
    %lowest price value for the year to ensure idle
    
    if counter == 0
            best_lambda_opt = min(price_data)-1;
    end
        
    % once lambda_opt has been found, take the unadjusted profit from that run
    [fake_prof, boost] = DAC_foropt(best_lambda_opt, price_data, X, k, parameters);
    best_prof = fake_prof-boost;
end