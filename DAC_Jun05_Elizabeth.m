
clc;

% Load pricing data
fileln = load('RTP_NYC_2010_2019.mat');
RTP = fileln.RTP;

% Set parameters
X_hat = 1; % Maximum ton capture per one absorption cycle
beta_a_1 = 0.2/X_hat; % Piecewise linear approximation for absorption, first-order term
beta_a_2 = -0.2/X_hat; % Piecewise linear approximation for absorption, second-order term
beta_d_1 = 0.0; % Piecewise linear approximation for desorption, first-order term
beta_d_2 = 0.4/X_hat; % Piecewise linear approximation for desorption, second-order term


% Initialization of price input (first day)
startingday = 1001;
L = RTP(1:T, startingday);

% action for a cycle of absorption (1) & desorption (0) takes length(Action) time periods.
Action = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0]; 
lambda_bar = 43; %threshold price
T = 288; %number of time periods per one simulation

%binary variables
u = true; % absorption on-off binary variable
v = false; % desorption on-off binary variable
k = false; %initial conditions
idle = false; %idle on-off binary variable
z = false; % switch cycle counting binary variable

   
for t = 1:T
    
    if t==1 %initial conditions
         k(t) = true; 
         z(t) = true;
    else
    
    if L(t) >= lambda_bar %if price is high, idle system
       u(t) = false;
       v(t) = false;
       
            if u(t-1) == true || idle(t-1) == true
                k(t) = true;
            else
                k(t) = false;
            end
            
       idle(t) = true;
      
    else
        idle(t) = false;
        J = mod(t-1, length(Action)) + 1; % takes remainder of t-1 divided
                                          % by length(Action) and adds 1
                                          % looks like this:
                                          % 1,2,3,4,5,6,7,8,9,10,1,2,3,4...
        if Action(J) == 1
            u(t) = true;
            v(t) = false;
            k(t) = 1;
        else
            u(t) = false;
            v(t) = true;
            k(t) = 0;
        end
        
        if k(t) == true && k(t-1) == false
         z(t) = true;
        else
         z(t) = false;
        end
        
    end
    end
    
    
end

% Initialize variables for a, d, and X
a = zeros(T, 1); %absorption amount of DAC system at time period t
d = zeros(T, 1); %desorption amount of DAC system at time period t
X = zeros(T, 1); %state-of-saturation capacity of DAC system at time period t

% Calculate values for a, d, and X
for t = 1:T
    
    if t==1 %initial conditions
         a(t) = beta_a_1 + beta_a_2 * 0;
         d(t) = beta_d_1 + beta_d_2 * 0;
         X(t) = 0+ a(t) - d(t); 
    else
        
        if u(t)
            a(t) = beta_a_1 + beta_a_2 * X(t-1);
        end
    
        if v(t)
            d(t) = beta_d_1 + beta_d_2 * X(t-1);
        end
    
        X(t) = X(t-1)+ a(t) - d(t); 
    
    end
end

fprintf(' t   u(t)   v(t)    k      z       L          X          a          d\n');
for t = 1:T 
    fprintf('%2d   %4d   %4d   %2d   %4d  %8.2f   %8.2f   %8.2f   %8.2f\n', t, u(t), v(t), k(t), z(t), L(t), X(t), a(t), d(t)); % Formatting and column widths
end


%problems

%1. a, d, and X are not between 0 and 1
