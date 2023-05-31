## package loading
using JuMP
using MAT
using Printf
using Gurobi
using CSV
using DataFrames
# load other packages before this line

@printf("Packages loaded")

## data loading 
# High-level Settings

# for this testing, using the same price data with battery arbitrage code
Zone = "NYC" # price zone name

# read price
fileln = matopen(string("./Compare_deg_model/RTP_",Zone,"_2010_2019.mat"))
RTP = read(fileln, "RTP")
close(fileln)
# This section loaded the price data in NYC from 2010-2019
# pricing data is the lambda_t in the formulation

## parameters setting

# maximium available DAC capacity 
X_hat = 1 # maximum X_hat ton capture per one absorption cycle

T = 72 # time step for optimization, 288 means 5 min temporal resolution for a day, 
# default: T = 288

S = 13 # paying S dollar per each cycle for material/operation cost
# default: S = 13

pi_co2 = 1000 # selling co2 at pi_co2 per ton
# default: pi_co2 = 180

P_a_unit = 0.642/19 # absorption power consumption MWh per time period per capacity
P_d_unit = 0.097/6 # desorption power consumption MWh per time period per capacity
# default: P_a_unit = 0.642/19
# default: P_d_unit = 0.097/6

P_a = P_a_unit*X_hat # correction of absorption power consumption by the installed capacity
P_d = P_d_unit*X_hat # correction of desoption power consumption by the installed capacity

beta_a_1 = 0.2/X_hat # peicewise linear approximation for absorption, first-order term
# default: beta_a_1 = 0.2/X_hat
beta_a_2 = -0.2/X_hat # peicewise linear approximation for absorption, second-order term
# default: beta_a_2 = -0.2/X_hat    !!!! beta_a_2 shall always be nagative of beta_a_1
beta_d_1 = 0.0 # peicewise linear approximation for desorption, first-order term
# default: beta_d_1 = 0.0           !!!! beta_d_1 shall always be 0 unless otherwise specified
beta_d_2 = 0.4/X_hat # peicewise linear approximation for desorption, second-order term
# default: beta_d_2 = 0.4/X_hat

# initialization for system
X0 = 0 # starting with no CO2 in the cycle
k0 = 0 # starting with no counting of the cycle

# initialization of price input (first day)
L = RTP[1:T,1];

# set a sufficiently large number for application of constraining binary variable related real numbers
M = 1000000000

## initialize optimization model

# setup model
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "MIPGap", 0.15) # set optimization attributes as needed
# set_silent(model) # no outputs

## set up variables

# binary variables
# absoprtion on-off binary variable
@variable(model, u[1:T], Bin)
# desorption on-off binary variable
@variable(model, v[1:T], Bin)
# sign variable for cycle counting
@variable(model, k[1:T], Bin)
# switch cycle counting binary variable
@variable(model, z[1:T], Bin)

# continuous variables
# state-of-saturation capacity variable
@variable(model, X[1:T], lower_bound = 0)
# absorption rate variable
@variable(model, a[1:T], lower_bound = 0)
# desorption rate variable
@variable(model, d[1:T], lower_bound = 0)

## constraints

# absorption/desorption binary constraint
@constraint(model, ab_vs_de[t=1:T], u[t] + v[t] <= 1)

# absorption/desorption rates constraints: initialize 
@constraint(model, ab_1, a[1] <= beta_a_1 + beta_a_2*((X[1]+X0)/2))
@constraint(model, de_1, d[1] <= beta_d_1 + beta_d_2*((X[1]+X0)/2))

# this initialization constraint shall be updated for each loop


# state-of-saturation capacity evolution constraint
@constraint(model, SOS_1, X[1] == a[1] - d[1] )
@constraint(model, SOS[t=2:T], X[t] - X[t-1] == a[t] - d[t])

# state-of-saturation capacity constraint
@constraint(model, SOS_U[t=1:T], X[t] <= X_hat)

# absorption/desorption rates constraints: general
@constraint(model, ab_t[t=2:T], a[t] <= beta_a_1 + beta_a_2*((X[t]+X[t-1])/2)) # absoprtion
@constraint(model, de_t[t=2:T], d[t] <= beta_d_1 + beta_d_2*((X[t]+X[t-1])/2)) # desorption

# absoprtion/desoroption rates constraints: bounded by binary variables
@constraint(model, ab_bi[t=1:T], a[t] <= u[t]*M)
@constraint(model, de_bi[t=1:T], d[t] <= v[t]*M)

# switch cycle constraint: initialize
@constraint(model, cycle_up_1, (-M)*(1-k[1]) <= k0 + (u[1] - v[1]) - 0.5)
@constraint(model, cycle_low_1, M*k[1] >= k0 + (u[1] - v[1]) - 0.5)
# switch cycle constraint: general
@constraint(model, cycle_up_t[t=2:T], (-M)*(1-k[t]) <= k[t-1] + (u[t] - v[t]) - 0.5)
@constraint(model, cycle_low_t[t=2:T], M*k[t] >= k[t-1] + (u[t] - v[t]) - 0.5)
# switch cycle constraint: identify the cycles
@constraint(model, cycle_check[t=2:T], z[t] >= k[t]-k[t-1])

# alternative cycle switching constraints
# @constraint(model, cycle_up[t=2:T], z[t] >= u[t]-u[t-1])
# @constraint(model, cycle_down[t=2:T], k[t] >= v[t]-v[t-1])
# !!! when using alternative switching constraints, objective function and set_objective_coefficient line shall be altered as well

# maximize revenue plus degradation value
@objective(model, Max, sum(pi_co2*d - L.*(P_a*u + P_d*v) - S*z))
# @objective(model, Max, sum(pi_co2*d - L.*(P_a*u + P_d*v) - S*z)) # this will triger very different results

# alternative form of objective function
# @objective(model, Max, sum(pi_co2*d - L.*(P_a*u + P_d*v) - S*(z+k)/2))


## initialize model for optimization
N_total = 3650; # total aviable number of simulations/days to run
N_sim = 5; # number of days, maximum 3650 for 10 years

# build up empty sets to store variables
R_s = zeros(1, N_sim) # profit (objective value) for each day of optimization
d_sim = zeros(1, N_sim) # captured CO2 for each round of optimization
L_S = zeros(T, N_sim) # get the price data as well
X_S = zeros(T, N_sim) # state_of_saturation for each day with T time steps
a_S = zeros(T, N_sim) # absoroption amount for each day with T time steps
d_S = zeros(T, N_sim) # desoroption amount for each day with T time steps
u_S = zeros(T, N_sim) # absoroption on-off status for each day with T time steps
v_S = zeros(T, N_sim) # desoroption on-off status for each day with T time steps
k_S = zeros(T, N_sim) # sign variable for each day with T time steps
z_S = zeros(T, N_sim) # cycle counting each day with T time steps

@printf("Optimization starts...\n")
@time begin
for n = (1):(N_sim)

    # update prices
    local L = RTP[1:T,n] # update the pricing data for each simulation to run 
    # update initial conditions
    # the initial condition for the next simulation takes last simulation's last value as initial input
    if n > 2
        local X0 = X_S[T, n-1]
        local k0 = k_S[T, n-1]
    end
    

    for t = 1:T
        # update the coefficient in the objective function
        set_objective_coefficient(model, u[t], -L[t]*P_a)
        set_objective_coefficient(model, v[t], -L[t]*P_d)
    end 

    optimize!(model)

    global R_s[n] = objective_value(model) # objective_value(model);
    global d_sim[n] = sum(value.(d)) # total captured CO2 per simulation

    global X_S[:,n] = value.(X) # state_of_saturation
    global a_S[:,n] = value.(a) # absorption amount
    global d_S[:,n] = value.(d) # desoroption amount
    global u_S[:,n] = value.(u) # absoprtion binary
    global v_S[:,n] = value.(v) # desoroption binary
    global k_S[:,n] = value.(k) # sign variable binary
    global z_S[:,n] = value.(z) # cycle counting binary

    # add the price data
    global L_S[:,n] = L 

    # termination_status(model)
    @printf("Finished Day %d, Cum Profit %d, Cum Tonnage %d, OptStatus: %s \n", n, sum(R_s), sum(d_S), termination_status(model))

end

end

results_sum = DataFrame(u = vec(u_S), v = vec(v_S), k = vec(k_S), z = vec(z_S), L = vec(L_S), X = vec(X_S), a = vec(a_S), d = vec(d_S)) # vectorize the results by columns and store 
show(results_sum)

results_obj = DataFrame(R = vec(R_s), D = vec(d_sim))


CSV.write("DAC_test_results.csv", results_sum) # write the results in to .csv file and save
CSV.write("DAC_revenue_results.csv", results_obj)