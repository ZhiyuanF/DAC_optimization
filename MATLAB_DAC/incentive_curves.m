clc
clearvars

% Set up length of simulation
steps_per_day = 288;
days = 365;

duration = steps_per_day * days;
startingday = 1;
start_step = steps_per_day * (startingday - 1) + 1;

% Define pi_co2 values
pi_co2_values = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500];

% Initialize results array (pi_co2 values, profit, net_co2)
results_array = zeros(length(pi_co2_values), 3);

%User input for setup
disp('Note: change the input pricing file name before running');
user_input1 = input('Enter "NY" for New York or "CA" for California to load emission data: ', 's');

% User input for adjusted or unadjusted CO2 price
adjust_choice = input('Do you want to use adjusted CO2 prices? (yes/no): ', 's');

tic

for idx = 1:length(pi_co2_values)
       pi_co2 = pi_co2_values(idx);
       disp(['Starting iteration for pi_co2 = ', num2str(pi_co2)]);  % Debugging statement

       %Parameters
        X_hat = 1;
        S = 115.6*X_hat;
        P_a_unit = 0.0642;
        P_d_unit = 0.02425;
        P_a = P_a_unit*X_hat;
        P_d = P_d_unit*X_hat;
        beta_a_1 = 0.2*X_hat;
        beta_a_2 = -0.2;
        beta_d_1 = 0.0*X_hat;
        beta_d_2 = 0.4;

        increment = 288;
        look_ahead = 288;

        parameters = [X_hat, S, pi_co2, P_a_unit, P_d_unit, P_a, P_d, beta_a_1, beta_a_2, beta_d_1, beta_d_2];

        %% Choose Increment and Look Ahead

        %Each “chunk” represents a portion of the data that is processed separately and has
        %its own optimized lambda. number of lambdas = number of chunks

        %Ex. 365 days, new lambda every 12 hrs, 24 hr look ahead:
        %   increment = 144
        %   look_ahead = 288
        %   number of lambdas (aka. number of chunks) = 365*2 = 730

        X = 0;
        k = 0;
        %% User Input Dependant Emission Loading

        %New York vs. California User Input for Emission Data

        
        if strcmpi(user_input1, 'NY')

            % Load emissions data from CSV
            myISO = readtable('NYISO_full_filtered_emission.csv');
            co2_data = transpose(myISO{1:duration, 11}); % Extract CO2 emissions from the 11th column  

        elseif strcmpi(user_input1, 'CA')

            % Load emissions data from CSV
            myISO = readtable('CAISO_full_filtered_emission.csv');
            co2_data = transpose(myISO{2:duration+1, 11});  % Extract CO2 emissions from the 11th column  

        else
            error('Invalid input. Please enter "NY" or "CA".');
        end
        
        %% Manual Price Loading

        % load the pricing data from the input .mat file (adjusted or unadjusted in
        % the title)
        input_mat_file = 'expanded_RAW_CAISO.mat';
        fileln = load(input_mat_file);
        myRTP = fileln.RTP;
        price_data = myRTP(start_step:start_step+duration-1); % Load prices for the specified time chunk

        % Initialize a tag for filenames
        adjust_tag = 'RAW';  % Default to 'RAW'

        % Check user choice and proceed accordingly
        if strcmpi(adjust_choice, 'yes')

            % Calculate adjusted CO2 price
            adjusted_price_data = price_data + pi_co2 * co2_data;

            % Update adjust_tag to include pi_co2 value
            adjust_tag = ['pi_' num2str(pi_co2)];

            % Save adjusted CO2 price into a .mat file
            adjusted_mat_filename = [adjust_tag '_CAISO.mat'];
            save(adjusted_mat_filename, 'adjusted_price_data');

            % Use adjusted_price_data for the rest of the code
            price_data = adjusted_price_data;

            myRTP = price_data;

        else
            % Use raw data
            myRTP = price_data;
        end

        %% Writing Final Test Sequence to CSV

        % Run LambdaOptimizer, get the optimal lambdas and profits for each chunk for
        % the whole data set
        [lambda_opts, profits] = LambdaOptimizer(price_data, increment, look_ahead, X, k, parameters);

        % create a cell with
        % # of rows = # of price data points
        % # of columns = 11
        final_results = cell(numel(price_data), 11);
       

        % loop through each chunk
        for chunk = 1:numel(lambda_opts)

            % length of the chunk is increment, except for at the very end of the
            % data set, where min rounds down
            chunk_len = min([increment, numel(price_data) - (chunk-1) * increment]);

            % column 1 is the length of the chunk
            %   every item of that column is = chunk number
            % column 2 is the length of the chunk
            %   every item of that column is the optimal lambda for that chunk
            c1 = cell(chunk_len,1);
            c1(:) = {chunk};
            c2 = cell(chunk_len,1);
            c2(:) = {lambda_opts{chunk}};

            % take the # of chunks that you're on, multiply by increment and then
            % add 1 to get the start point
            % take the # of chunks that you're on, multiply by increment and then
            % add the chunk length to get the end point
            final_results((chunk-1) * increment +1 : (chunk-1) * increment + chunk_len, 1) = c1(:);
            final_results((chunk-1) * increment +1 : (chunk-1) * increment + chunk_len, 2) = c2(:);

            % retrieve right start and end point for the chunk in price_data
            small_price_data = price_data(1, (chunk-1) * increment +1 : (chunk-1) * increment + chunk_len);

            % columns 3-11 come from DAC_fordata function that runs on whatever the
            % current optimal lambda value is for that chunk

            % note: line 133 is what takes the longest, it's not the optimization
            % it's the data writing
            final_results((chunk-1) * increment +1 : (chunk-1) * increment + chunk_len, 3:11) = DAC_fordata(lambda_opts{chunk}, small_price_data, X, k, parameters);

            % update X and k at the beginning of the next chunk
            [X, k] = DAC_foriteration(lambda_opts{chunk}, small_price_data, X, k, parameters);
        end

        disp(sum([profits{:}]))

        toc

        %% Writing Day-by-day results to CSV

        % total_profits_per_day is a cell array that stores the total profit for
        % each 288-step window, or hoever long a day is defined to be

        total_profits_per_day = cell(365, 4);

        % Initialize arrays to store the cumulative CO2 captured (both gross and net)
        cumulative_co2_captured_gross = zeros(days, 1);
        cumulative_co2_captured_net = zeros(days, 1);

        % loop through each day REGARDLESS OF CHUNK SIZE
        for day = 1:days
            start_row = (day - 1) * steps_per_day + 1;
            end_row = start_row + steps_per_day - 1;

            %store profit, gross CO2, emitted CO2, and net CO2
            profits_per_day = cell2mat(final_results(start_row:end_row, 11));
            desorb_co2_day = sum(cell2mat(final_results(start_row:end_row, 10)));
            co2_emitted_day = P_a * dot(cell2mat(final_results(start_row:end_row, 3)), co2_data(start_row:end_row)') + P_d * dot(cell2mat(final_results(start_row:end_row, 4)), co2_data(start_row:end_row)');
            net_co2_day = desorb_co2_day - co2_emitted_day;

            %writing to total_profits_per_day
            total_profits_per_day{day, 1} = startingday + day - 1;
            total_profits_per_day{day, 2} = sum(profits_per_day);
            total_profits_per_day{day, 3} = desorb_co2_day;
            total_profits_per_day{day, 4} = net_co2_day;
        end
        
        % calculate the sum of profits and store it in results array
        profits_sum = sum([profits{:}]);
        results_array(idx, 2) = profits_sum;

        % calculate the net CO2 and store it in results array
        net_co2 = total_profits_per_day(:, 4);
        net_co2_numeric = cell2mat(net_co2); % convert cell to numeric array
        results_array(idx, 3) = sum(net_co2_numeric);
end

% add pi_co2 values in first column
results_array(:, 1) = pi_co2_values;

% write to csv
writematrix(results_array, 'sensitivity_Extended_adjusted_caiso.csv');

%end
toc
