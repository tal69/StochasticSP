import os

import numpy as np
import Single_SP
import subprocess
import time
import datetime


"""  Prepare and run all the three models version with closet SP only """
"""  June 26, 2022  by Tal Raviv """

""" User input parameters """
replications_range = range(1,10)  # number of random replications for each settings combination of m,n,r, and p
C = np.array( [30,50, 80])   # Possible SP capacities
C_num = len(C)
alpha_range = [50, 200] # range of the weights of rejectionsin the objective function
p_range = [0.3, 0.5]   # range of periodic pikcup probability
m_range = [50, 100]      # range of the number of candidate SP locations
n_range = [200,400]   # range of the number of demand point in the models
r_range = [0.1, 0.15]  # Service radius range (since the locations are always generated on a unit square this value shudld decrease when m and n are increased
safety_margin_range = [1]  # set to empty list if there is no need to solve the deterministic model
scenario_num = 30 # Number of scenarios to generate and use in the stochastic programing model - set to 0 if the scenario model should not be run
time_limit = 900  # time limit to be passed to the solver (seconds)
do_our = True  # set False to skip our model


"""  The values of mu, H, Df, and Fd are randomly generated based on the above parameters """

file_name  = "Results_SP_v5.csv"  # Don't change this, because the OPL model assumes this name

""" definition for generating Lambda_k """
rho = [0, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.25, 1.5,  2]
Lambda_points = len(rho)
f = open(file_name, "a")

now = datetime.datetime.now()

f.write(f"\n rep,n,m,p, r, alpha, xs,Total demand")
if do_our:
    f.write(f",cplex status, ObjValue, LB, Cplex Time, SP Cost, model rej.")
    for s in range(C_num):
        f.write(f", SPs ({C[s]})")

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
        for m in m_range:
            sp_locations = np.random.rand(m,2)  # SP locations

            h = np.int_(5+np.dot(np.random.normal(p*5,p,m).reshape(m,1), 1+C.reshape(1,C_num))) # SP costs
            for n in n_range:
                """   Generate demand for each each demand point   """
                mu = np.random.randint(50,int(2000*m/n), size=n) / 100  # the demand of each demand point
                for r in r_range:
                    polar = np.random.rand(n,2)
                    polar[:,0] *= 2*np.pi
                    polar[:,1] *= r
                    demand_points = np.zeros((n,2))
                    for d in range(n):
                        f = d % m
                        demand_points[d,0] = sp_locations[f,0] + r * np.cos(polar[d,0]) * polar[d,1]
                        demand_points[d,1] = sp_locations[f,1] + r * np.sin(polar[d,0]) * polar[d,1]

                    # Count the number of actual x variables in the model
                    number_of_xs = sum([ np.linalg.norm((demand_points[d] - sp_locations[f]), ord=2) < r for d in range(n) for f in range(m) ])

                    for alpha in alpha_range:
                        f = open(file_name, "a")
                        f.write(f"{rep}, {n},{m},{p}, {r}, {alpha}, {number_of_xs}, {np.sum(mu)}")
                        f.close()

                        if do_our:

                            f = open('SP_loc_cap_v5.dat', "w")
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
                                subprocess.run(["oplrun", "sp_loc_cap_v5.mod", "SP_loc_cap_v5.dat"])
                            except:
                                print("Can't run OPL model (sp_loc_cap_v5)")
                                exit(1)

                        for safety_margin_idx in range(len(safety_margin_range)):
                            safety_margin = safety_margin_range[safety_margin_idx]
                            f = open('SP_cover_v5.dat', "w")
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
                                subprocess.run(["oplrun", "sp_cover_v5.mod", "SP_cover_v5.dat"])
                            except:
                                print("Can't run OPL model (sp_cover_v5)")
                                exit(1)

                            if not os.path.exists("out_cover.txt"):
                                print("Panic: oplrun failed to produce output file")
                                exit(1)

                            " Read solution of the cover model"
                            f = open("out_cover.txt","r")
                            ez = [ float(s) for s in f.readlines()]
                            sp_demand = ez[:m]
                            sp_type = np.array([int(i) for i in ez[m:]])
                            f.close()

                            " estimate number of rejections in solution of the covering model "
                            Rej = 0
                            for i in range(m):
                                if sp_type[i]:
                                    Rej += np.interp(sp_demand[i]/safety_margin, Lambda[(sp_type[i]-1),:], R[(sp_type[i]-1),:])

                            " CHECK HERE - maybe we don't count uncovered SPs  "

                            f = open(file_name, "a")
                            f.write(f",{Rej}")
                            f.close()

                        if scenario_num > 0:
                            f = open('SP_scenarios_v5.dat', "w")
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
                                subprocess.run(["oplrun", "sp_scenarios_v5.mod", "SP_scenarios_v5.dat"])
                            except:
                                print("Can't run OPL model (SP_scenarios_v5)")
                                exit(1)

                            if not os.path.exists("out_scenarios.txt"):
                                print("Panic: oplrun failed to produce output file (out_scenarios.txt)")
                                exit(1)

                            " Read solution of the cover model"
                            f = open("out_scenarios.txt","r")
                            sp_type  = [ int(s) for s in f.readlines()]
                            f.close()

                            """ find closest open SP of each demand point """
                            closest_sp = -np.ones(n, int)
                            closest_sp_dist = np.ones(n)   # relevant distance are always less than r
                            for f in range(m):
                                if sp_type[f] != 0:
                                    for d in range (n):
                                        dist = np.linalg.norm((demand_points[d] - sp_locations[f]), ord=2)
                                        if dist <= closest_sp_dist[d]:
                                            closest_sp_dist[d] = dist
                                            closest_sp[d] = f


                            """ calculate demand at each SP """
                            Rej = 0
                            sp_demand = np.zeros(m)
                            for d in range(n):
                                if closest_sp_dist[d] > r:
                                    print(f"Panic: closest SP of  {d} is {closest_sp[d]} at distance {closest_sp_dist[d]}")
                                    exit(1)
                                sp_demand[closest_sp[d]] += mu[d]

                            " estimate number of rejections in the solution of the scenario base model model "

                            for i in range(m):
                                if sp_type[i]:
                                    Rej += np.interp(sp_demand[i], Lambda[(sp_type[i]-1),:], R[(sp_type[i]-1),:])

                            f = open(file_name, "a")
                            f.write(f",{Rej}")
                            f.close()


                        f = open(file_name, "a")
                        f.write("\n")
                        f.close()
