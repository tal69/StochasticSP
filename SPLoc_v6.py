import os

import numpy as np
import Single_SP
import subprocess
import time


"""  Prepare and run all the three models version with closet SP only """
""" Same as v5 but with locations on a gird and with post calculation of the objective function """
"""  July 2, 2022  by Tal Raviv """

""" User input parameters """
replications_range = range(10)  # number of random replications for each settings combination of m,n,r, and p
C = np.array( [30,60,90])   # Possible SP capacities
C_num = len(C)
alpha_range = [5,20] # range of the weights of rejectionsin the objective function
p_range = [0.5]   # range of periodic pikcup probability
demand_points_sqr = 21
demand_points_dist = 200
sp_sqr = 10
sp_dist = 400
setup_cost_fix = 10
setup_cost_locker = 2/3  # The cost of a locker is sampled from (setup_cost_fix + setup_cost_locker*C ) * normal(1,1/6)

n = demand_points_sqr ** 2
m = sp_sqr ** 2

r_range = [401, 601]  # Service radius range (since the locations are always generated on a unit square this value shudld decrease when m and n are increased
safety_margin_range = [1]  # set to empty list if there is no need to solve the deterministic model
scenario_num = 50 # Number of scenarios to generate and use in the stochastic programing model - set to 0 if the scenario model should not be run
time_limit = 90  # time limit to be passed to the solver (seconds)
do_our = True  # set False to skip our model


def calc_objective_value(sol_file_name, Lambda, R, sp_locations, demand_points, mu,r):

    if not os.path.exists(sol_file_name):
        print(f"Panic: oplrun failed to produce output file ({sol_file_name})")
        exit(1)

    m = sp_locations.shape[0]
    n = demand_points.shape[0]

    " Read solution of the cover model"
    ff = open(sol_file_name, "r")
    ez = [float(s) for s in ff.readlines()]
    ff.close()
    sp_type = np.array([int(i) for i in ez[:m]])
    sp_demand_model = ez[m:]

    """  The file may contain also demand info but we don't use it now """


    """ find closest open SP of each demand point """
    closest_sp = -np.ones(n, int)
    closest_sp_dist = np.ones(n)*r*2   # relevant distance are always less than r
    for sp in range(m):
        if sp_type[sp] != 0:
            for dp in range (n):
                dist = np.linalg.norm((demand_points[dp] - sp_locations[sp]), ord=2)
                if dist <= closest_sp_dist[dp]:
                    closest_sp_dist[dp] = dist
                    closest_sp[dp] = sp


    """ calculate demand at each SP """
    rejections = 0
    sp_demand = np.zeros(m)
    for dp in range(n):
        if closest_sp_dist[dp] > r:
            print(f"Panic: closest SP of  {dp} is {closest_sp[dp]} at distance {closest_sp_dist[dp]}")
            exit(1)
        sp_demand[closest_sp[dp]] += mu[dp]

    " estimate number of rejections in the solution of the scenario base model model "

    for i in range(m):
        if sp_type[i]:
            rejections += np.interp(sp_demand[i], Lambda[(sp_type[i]-1),:], R[(sp_type[i]-1),:])

    return rejections


"""  The values of mu, H, Df, and Fd are randomly generated based on the above parameters """

file_name = "Results_SP_v6.csv"  # Don't change this, because the OPL model assumes this name

""" definition for generating Lambda_k """
rho = [0, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.25, 1.5,  2]
Lambda_points = len(rho)
f = open(file_name, "a")

f.write(f"\n rep,n,m,p, r, alpha, xs,Total demand")
if do_our:
    f.write(f",cplex status, ObjValue, LB, Cplex Time, SP Cost, model rej.")
    for s in range(C_num):
        f.write(f", SPs ({C[s]})")
    f.write(", Rejections Calc")

for safety_margin in safety_margin_range:
    f.write(f",safety_margin, cplex status, ObjValue, LB, Cplex Time, SP Cost, Rejections Model")
    for s in range(C_num):
        f.write(f", SPs ({C[s]})")
    f.write(", Rejections Calc")

if scenario_num:
    f.write(f",scenarios_num, cplex status, ObjValue, LB, Cplex Time, SP Cost, model rej.")
    for s in range(C_num):
        f.write(f", SPs ({C[s]})")
    f.write(", Rejections Calc")

f.write("\n")
f.close()

np.set_printoptions(threshold=np.inf)  # cancel numpy truncation when printing large arrays

demand_points = np.zeros((n,2))
sp_locations = np.zeros((m,2))
d = 0
for x in range(demand_points_dist // 2, demand_points_sqr*demand_points_dist, demand_points_dist):
    for y in range(demand_points_dist // 2, demand_points_sqr*demand_points_dist, demand_points_dist):
        demand_points[d,:]  = (x,y)
        d += 1

for p in p_range:# time-to-pickup distribution

    Lambda = np.zeros((C_num, Lambda_points))
    R = np.zeros((C_num, Lambda_points))
    print(f"Calculating Lambda for p={p}")
    start_time = time.time()
    for s in range(C_num):
        for k in range(Lambda_points):
            Lambda[s,k] = rho[k] * C[s] * p
            R[s,k] = Single_SP.markov_chain_sp(C[s], Lambda[s,k], p)[0]

    for rep in replications_range:
        np.random.seed(rep)

        s = 0
        for x in range((sp_dist+demand_points_dist) // 2, sp_sqr*sp_dist, sp_dist):
            for y in range(sp_dist // 2, sp_sqr*sp_dist, sp_dist):
                sp_locations[s,:] = (x + np.random.normal(0, demand_points_dist/12),y + np.random.normal(0, demand_points_dist/12))
                s += 1

        #h = np.int_(np.dot(np.random.normal(1,1/6,m).reshape(m,1), sprint("Status:", LpStatus[prob.status])hape(1,C_num))) # SP costs
        h = np.int_(np.dot(np.random.normal(1,1/6,m).reshape(m,1), setup_cost_fix+setup_cost_locker*C.reshape(1,C_num))) # SP costs

        """   Generate demand for each each demand point   """
        mu = np.random.randint(50,1000, size=n) / 100  # the demand of each demand point
        for r in r_range:

            # Count the number of actual x variables in the model
            number_of_xs = sum([ np.linalg.norm((demand_points[d] - sp_locations[f]), ord=2) < r for d in range(n) for f in range(m) ])

            for alpha in alpha_range:
                f = open(file_name, "a")
                f.write(f"{rep}, {n},{m},{p}, {r}, {alpha}, {number_of_xs}, {np.sum(mu)}")
                f.close()

                if do_our:
                    f = open('SP_loc_cap_v6.dat', "w")
                    f.write(f"time_limit = {time_limit};\n")
                    f.write(f"alpha = {alpha};\n")
                    f.write(f"r = {r};\n")
                    f.write(f"n = {n};\n")
                    f.write(f"m = {m};\n")
                    f.write(f"C_num = {C_num};\n")
                    f.write(f"h = {h};\n")
                    f.write(f"mu = {mu};\n")
                    f.write(f"Lambda_points = {Lambda_points};\n")
                    f.write(f"Lambda = {Lambda};\n")
                    f.write(f"R = {R};\n")
                    f.write(f"SP_xy = {sp_locations};\n")
                    f.write(f"DP_xy = {demand_points};\n")
                    f.close()

                    try:
                        os.remove("out_storage.txt")   # delete so we know if oplrun produced new output
                    except:
                        print("Can't delete out_storage.txt. Maybe because of Dropbox or because it is the first iteration and the file was not created yet. Let us hope for the best")

                    try:
                        subprocess.run(["oplrun", "sp_loc_cap_v6.mod", "SP_loc_cap_v6.dat"])
                    except:
                        print("Can't run OPL model (sp_loc_cap_v6)")
                        exit(1)

                    f = open(file_name, "a")
                    f.write(f",{calc_objective_value('out_our.txt', Lambda, R, sp_locations, demand_points, mu,r)}")
                    f.close()

                for safety_margin_idx in range(len(safety_margin_range)):
                    safety_margin = safety_margin_range[safety_margin_idx]
                    f = open('SP_cover_v6.dat', "w")
                    f.write(f"time_limit = {time_limit};\n")
                    f.write(f"alpha = {alpha};\n")
                    f.write(f"r = {r};\n")
                    f.write(f"n = {n};\n")
                    f.write(f"m = {m};\n")
                    f.write(f"C_num = {C_num};\n")
                    f.write(f"C = {C};\n")
                    f.write(f"h = {h};\n")
                    f.write(f"mu = {mu};\n")
                    f.write(f"p = {p};\n")
                    f.write(f"safety_margin = {safety_margin};\n")
                    f.write(f"SP_xy = {sp_locations};\n")
                    f.write(f"DP_xy = {demand_points};\n")
                    f.close()

                    try:
                        os.remove("out_cover.txt")   # delete so we know if oplrun produced new output
                    except:
                        print("Can't delete out_cover.txt (maybe because of Dropbox or because it is the first iteration and out_cover.txt was not created yet) let us hope for the best")

                    try:
                        subprocess.run(["oplrun", "sp_cover_v6.mod", "SP_cover_v6.dat"])
                    except:
                        print("Can't run OPL model (sp_cover_v6)")
                        exit(1)

                    f = open(file_name, "a")
                    f.write(f",{calc_objective_value('out_cover.txt', Lambda, R, sp_locations, demand_points, mu,r)}")
                    f.close()

                if scenario_num > 0:
                    f = open('SP_scenarios_v6.dat', "w")
                    f.write(f"time_limit = {time_limit};\n")
                    f.write(f"alpha = {alpha};\n")
                    f.write(f"r = {r};\n")
                    f.write(f"n = {n};\n")
                    f.write(f"m = {m};\n")
                    f.write(f"C_num = {C_num};\n")
                    f.write(f"p = {p};\n")
                    f.write(f"C = {C};\n")
                    f.write(f"h = {h};\n")
                    f.write(f"scenario_num={scenario_num};\n")
                    f.write(f"mu = {np.random.poisson(mu, size = (scenario_num,n)).reshape(n,scenario_num)};\n")
                    f.write(f"SP_xy = {sp_locations};\n")
                    f.write(f"DP_xy = {demand_points};\n")
                    f.close()

                    try:
                        os.remove("out_scenarios.txt")   # delete so we know if oplrun produced new output
                    except:
                        print("Can't delete out_cover.txt (maybe because of Dropbox) let us hope for the best")

                    try:
                        subprocess.run(["oplrun", "sp_scenarios_v6.mod", "SP_scenarios_v6.dat"])
                    except:
                        print("Can't run OPL model (SP_scenarios_v6)")
                        exit(1)

                    if not os.path.exists("out_scenarios.txt"):
                        print("Panic: oplrun failed to produce output file (out_scenarios.txt)")
                        exit(1)

                    f = open(file_name, "a")
                    f.write(f",{calc_objective_value('out_scenarios.txt', Lambda, R, sp_locations, demand_points, mu,r)}")
                    f.close()

                f = open(file_name, "a")
                f.write("\n")
                f.close()
