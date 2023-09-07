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
# fileln = matopen(string("./Compare_deg_model/RTP_",Zone,"_2010_2019.mat"))
fileln = matopen(string("./Compare_deg_model/RTP_",Zone,"_2010_2019.mat"))
RTP = read(fileln, "RTP")
close(fileln)
# This section loaded the price data in NYC from 2010-2019
# pricing data is the lambda_t in the formulation

RTPv = vec(RTP) # convert matrix to vector
n_N = 288 # there are 288 steps per day

## parameters setting

# maximium available DAC capacity 
X_hat = 10 # maximum X_hat ton capture per one absorption cycle

T = 24 # optimization horizon
H = 12 # time period record optimization results
nH = Int(288/H) # number of optimization periods each day

S = 115.6*X_hat # paying S dollar per each cycle for material/operation cost
# default: S = 115.6

pi_co2 = 200 # selling co2 at pi_co2 per ton
# default: pi_co2 = 180

P_a_unit = 0.642/10 # absorption power consumption MWh per time period per capacity
P_d_unit = 0.097/4 # desorption power consumption MWh per time period per capacity
# default: P_a_unit = 0.642/19
# default: P_d_unit = 0.097/6

P_a = P_a_unit*X_hat # correction of absorption power consumption by the installed capacity
P_d = P_d_unit*X_hat # correction of desoption power consumption by the installed capacity

beta_a_1 = 0.2*X_hat # peicewise linear approximation for absorption, first-order term
# default: beta_a_1 = 0.2*X_hat
beta_a_2 = -0.2 # peicewise linear approximation for absorption, second-order term
# default: beta_a_2 = -0.2  
beta_d_1 = 0.0*X_hat # peicewise linear approximation for desorption, first-order term
# default: beta_d_1 = 0.0*X_hat           
beta_d_2 = 0.4 # peicewise linear approximation for desorption, second-order term
# default: beta_d_2 = 0.4

# initialization for system
X0 = 0 # starting with no CO2 in the cycle
k0 = 0 # starting with no counting of the cycle

# initialization of price input (first day)
L = RTP[1:T,1];

# set a sufficiently large number for application of constraining binary variable related real numbers
M = 10

## initialize optimization model

# setup model
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "MIPGap", 0.05) # set optimization attributes as needed
set_silent(model) # no outputs

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


# this initialization constraint shall be updated for each loop


# state-of-saturation capacity evolution constraint
@constraint(model, SOS_1, X[1] - X0 == a[1] - d[1] )
@constraint(model, SOS[t=2:T], X[t] - X[t-1] == a[t] - d[t])

# state-of-saturation capacity constraint
@constraint(model, SOS_U[t=1:T], X[t] <= X_hat)

# absorption/desorption rates constraints: initialize 
@constraint(model, ab_1, a[1] <= beta_a_1 + beta_a_2*X0)
@constraint(model, de_1, d[1] <= beta_d_1 + beta_d_2*X0)
# absorption/desorption rates constraints: general
@constraint(model, ab_t[t=2:T], a[t] <= beta_a_1 + beta_a_2*X[t-1]) # absoprtion
@constraint(model, de_t[t=2:T], d[t] <= beta_d_1 + beta_d_2*X[t-1]) # desorption

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
@constraint(model, cycle_check_1, z[1] >= k[1]-k0)
@constraint(model, cycle_check[t=2:T], z[t] >= k[t]-k[t-1])

# alternative cycle switching constraints
# @constraint(model, cycle_up[t=2:T], z[t] >= u[t]-u[t-1])
# @constraint(model, cycle_down[t=2:T], k[t] >= v[t]-v[t-1])
# !!! when using alternative switching constraints, objective function and set_objective_coefficient line shall be altered as well

# maximize revenue plus degradation value
# @objective(model, Min, -sum(pi_co2*d - L.*(P_a*u + P_d*v) - S*z))
@objective(model, Max, 100+sum(pi_co2*d - L.*(P_a*u + P_d*v) - S*z))

# alternative form of objective function
# @objective(model, Max, sum(pi_co2*d - L.*(P_a*u + P_d*v) - S*(z+k)/2))


## initialize model for optimization
N_total = 3650; # total available number of simulations/days to run
N_start = 1119; # number of starting day
N_sim = 5; # number of days, maximum 3650 for 10 years

# build up empty sets to store variables
R_s = zeros(1, N_sim) # profit (objective value) for each day of optimization
d_sim = zeros(1, N_sim) # captured CO2 for each day of optimization
L_S = zeros(n_N, N_sim) # get the price data as well
X_S = zeros(n_N, N_sim) # state_of_saturation for each day with T time steps
a_S = zeros(n_N, N_sim) # absoroption amount for each day with T time steps
d_S = zeros(n_N, N_sim) # desoroption amount for each day with T time steps
u_S = zeros(n_N, N_sim) # absoroption on-off status for each day with T time steps
v_S = zeros(n_N, N_sim) # desoroption on-off status for each day with T time steps
k_S = zeros(n_N, N_sim) # sign variable for each day with T time steps
z_S = zeros(n_N, N_sim) # cycle counting each day with T time steps

@printf("Optimization starts...\n")
@time begin
for n = 1:N_sim # go through each day

for t = 1:nH# go through each step
 
    # update initial conditions
    # the initial condition for the next simulation takes last simulation's last value as initial input
    if n > 1 && t == 1
        local X0 = X_S[n_N, n-1]
        local k0 = k_S[n_N, n-1]

        set_normalized_rhs(ab_1, beta_a_1 + beta_a_2*X0)
        set_normalized_rhs(de_1, beta_d_1 + beta_d_2*X0)
        set_normalized_rhs(SOS_1, X0)
        set_normalized_rhs(cycle_up_1, M+k0-0.5)
        set_normalized_rhs(cycle_low_1, k0-0.5)
        set_normalized_rhs(cycle_check_1, -k0)
    end

    # update constraint
    if t > 1
        local X0 = X_S[Int((t-1)*H), n]
        local k0 = k_S[Int((t-1)*H), n]

        set_normalized_rhs(ab_1, beta_a_1 + beta_a_2*X0)
        set_normalized_rhs(de_1, beta_d_1 + beta_d_2*X0)
        set_normalized_rhs(SOS_1, X0)
        set_normalized_rhs(cycle_up_1, M+k0-0.5)
        set_normalized_rhs(cycle_low_1, k0-0.5)
        set_normalized_rhs(cycle_check_1, -k0)
    end      

    # Update objective
    for tau = 1:T
        # update the coefficient in the objective function
        set_objective_coefficient(model, u[tau], -RTPv[Int((n+N_start-2)*n_N + tau + (t-1)*H)]*P_a)
        set_objective_coefficient(model, v[tau], -RTPv[Int((n+N_start-2)*n_N + tau + (t-1)*H)]*P_d)
    end 

    optimize!(model)

    @printf("Finished Period %d of Day %d, Obj: %.2f, OptStatus: %s \n", t, n, objective_value(model), termination_status(model))


    # global R_s[n] = objective_value(model) # objective_value(model);
    # global d_sim[n] = sum(value.(d)) # total captured CO2 per simulation

    local t2H = Int((t-1)*H)

    global X_S[t2H + 1:t2H + H,n] = value.(X[1:H]) # state_of_saturation
    global a_S[t2H + 1:t2H + H,n] = value.(a[1:H]) # absorption amount
    global d_S[t2H + 1:t2H + H,n] = value.(d[1:H]) # desoroption amount
    global u_S[t2H + 1:t2H + H,n] = value.(u[1:H]) # absoprtion binary
    global v_S[t2H + 1:t2H + H,n] = value.(v[1:H]) # desoroption binary
    global k_S[t2H + 1:t2H + H,n] = value.(k[1:H]) # sign variable binary
    global z_S[t2H + 1:t2H + H,n] = value.(z[1:H]) # cycle counting binary



end
    # update prices
    local L = RTP[:,n+N_start-1] # update the pricing data for each simulation to run 

    global R_s[n] = sum(pi_co2*d_S[:,n] - L.*(P_a*u_S[:,n] + P_d*v_S[:,n]) - S*z_S[:,n]) # objective_value(model);
    global d_sim[n] = sum(d_S[:,n]) # total captured CO2 per simulation
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