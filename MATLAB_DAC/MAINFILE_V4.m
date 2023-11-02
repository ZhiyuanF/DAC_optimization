clc
clearvars

%Set up length of simulation
steps_per_day = 288;
days = 365;

duration = steps_per_day * days;
startingday = 1;
start_step = steps_per_day * (startingday - 1) + 1;

pi_co2 = 200; % selling co2 at pi_co2 per ton

%Parameters
X_hat = 23.95;
S = 213.36*X_hat;
P_a_unit = 0.000353;
P_d_unit = 0.001427;
P_a = P_a_unit*X_hat;
P_d = P_d_unit*X_hat;
beta_a_1 = 0.00099*X_hat;
beta_a_2 = 0;
beta_d_1 = 0;
beta_d_2 = 0.088;
increment = 3211;
look_ahead = 3211;

parameters = [X_hat, S, pi_co2, P_a_unit, P_d_unit, P_a, P_d, beta_a_1, beta_a_2, beta_d_1, beta_d_2];

%User input for setup
disp('Note: change the input pricing file name before running');
user_input1 = input('Enter "NY" for New York or "CA" for California to load emission data: ', 's');

tic

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
         
    % Load emissions data from CSVN
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

% extract the base filename without the extension from the loaded .mat file
[~, base_filename, ~] = fileparts(input_mat_file);

% if the base_filename starts with 'RTP_', remove it
if strncmp(base_filename, 'RTP_', 4)
    base_filename = base_filename(5:end);
end

% construct the output CSV filenames based on the input dataset
results_filename = ['MAT_' base_filename '_test.csv'];
total_profits_filename = ['MAT_' base_filename '_revenue.csv'];

myRTP = price_data; % Assign the transposed price data to myRTP

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

% filename for sequence printing
fileID = fopen(results_filename, 'w');

% column names to the sequence CSV file
fprintf(fileID, 'chunk,optimal lambda,u,v,k,z,L,X,a,d,profit\n');

% loop through the cell array and write each row as a CSV line
for row = 1:size(final_results, 1)
    fprintf(fileID, '%d,%f,%.0f,%.0f,%.0f,%.0f,%f,%f,%f,%f,%f\n', final_results{row,:});
end

% close the file
fclose(fileID);

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
    
    % store cumulative gross and net CO2
    cumulative_co2_captured_gross(day) = sum(cell2mat(final_results(1:end_row, 10))); 
    cumulative_co2_captured_net(day) = sum(cell2mat(final_results(1:end_row, 10))) - (P_a * dot(cell2mat(final_results(1:end_row, 3)), co2_data(1:end_row)') + P_d * dot(cell2mat(final_results(1:end_row, 4)), co2_data(1:end_row)')); 
end

% filename for total profits per day
fileID_total_profits = fopen(total_profits_filename, 'w');

% column names to the CSV file
fprintf(fileID_total_profits, 'day,profit,CO2_captured,net_CO2_day\n');

% loop through the cell array and write each row as a CSV line
for row = 1:size(total_profits_per_day, 1)
fprintf(fileID_total_profits, '%d,%f,%f,%f\n', total_profits_per_day{row, :});
end 

% close the file
fclose(fileID_total_profits);

%end
toc
