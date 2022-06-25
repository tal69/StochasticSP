import numpy as np
import Single_SP
import subprocess
import time
import datetime


""" User input paramteres - for the experiement in Section 5.2
    Warning this may run for up to 240 Hours!
"""

replications = 10  # number of random replications for each settings combination of m,n,r, and p
C = np.array( [20,40, 60, 80, 100])   # Possible SP capacities
C_num = len(C)
alpha = 10
p_range = [0.5, 0.7]
m_range = [300,500,1000]
n_range = [2000]
r_range = [0.05, 0.15]
safety_margin_range = []  # set to empty list if there is no need to solve the covering model

"""  The values of mu, H, Df, and Fd are randomly generetaed based on the above paramtered"""

file_name = "Results_SP.csv"  # Don't change this, beacuse the OPL model assumes this name

""" definition for genrating Lambda_k """
Lambda_max_ratio = 3
Lambda_min_ratio = 0.3
Lambda_points = 11

f = open(file_name, "a")

now = datetime.datetime.now()

#f.write(f"\nStart time: {now.strftime('%Y-%m-%d %H:%M:%S)')}\n")
f.write(f"\nrep,n,m,p, r, xs,Total demand, cplex status, ObjValue, LB, Cplex Time, SP Cost")
for s in range(C_num):
    f.write(f", SPs ({C[s]})")
f.write(", Rejections")
for safety_margin in safety_margin_range:
    f.write(f",safety_margin, cplex status, ObjValue, LB, Cplex Time")
    for s in range(C_num):
        f.write(f", SPs ({C[s]})")
    f.write(", Rejections")
f.write("\n")
f.close()


np.set_printoptions(threshold=np.inf)  # cancel numpy truncation when printing large arrays

for p in p_range:# time-to-pickup distribution

    Lambda = np.zeros((C_num, Lambda_points))
    R = np.zeros((C_num, Lambda_points))
    print(f"Calcuating Lambda for p={p}")
    start_time = time.time()
    for s in range(C_num):
        Lambda_min = Lambda_min_ratio *  p * C[s]
        Lambda_high = Lambda_max_ratio *  p * C[s]
        for k in range(1,Lambda_points):
            Lambda[s,k] = Lambda_min+(k-1)*(Lambda_high-Lambda_min)/(Lambda_points-2)
            R[s,k] = Single_SP.markov_chain_sp(C[s], Lambda[s,k], p)[0]
            """ Lambda[s,0] and R[s,0] were left 0 intentionally """
    f = open(file_name, "a")
    f.write(f"Preprocessing time for p={p}:  {time.time()-start_time}\n")
    f.close()

    for rep in range(replications):
        np.random.seed(rep)
        for m in m_range:
            sp_locations =  np.random.rand(m,2)
            " Generate SP costs "
            h = np.int_(5+np.dot(np.random.normal(alpha*p*0.5,alpha*p*0.1,m).reshape(m,1), 1+C.reshape(1,C_num)))
            for n in n_range:
                "  Generate demand for each each demaned point   "
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
                    Fd = {}
                    Df = {}

                    for d in range(n):
                        Fd[d] = set([])

                    for f in range(m):
                        Df[f] = set([])

                    number_of_xs = 0
                    for f in range(m):
                        for d in range(n):
                            if np.linalg.norm((demand_points[d] - sp_locations[f]), ord=2) < r:
                                Df[f].add(d)
                                Fd[d].add(f)
                                number_of_xs += 1

                    f = open('SP_loc_cap.dat', "w")
                    f.write(f"alpha = {alpha};\n")
                    f.write(f"n = {n};\n")
                    f.write(f"m = {m};\n")
                    f.write(f"C_num = {C_num};\n")
                    f.write(f"h = {h};\n")
                    f.write(f"mu = {mu};\n")
                    f.write(f"Lambda_points = {Lambda_points};\n")
                    f.write(f"Lambda = {Lambda};\n")
                    f.write(f"R = {R};\n")
                    f.write(f"Fd = {list(Fd.values())};\n")
                    f.write(f"Df = {list(Df.values())};\n")

                    """ For reporting only """
                    f.close()

                    f = open(file_name, "a")
                    f.write(f"{rep}, {n},{m},{p}, {r}, {number_of_xs}, {np.sum(mu)},")
                    f.close()

                    try:
                        subprocess.run(["oplrun", "sp_loc_cap.mod", "SP_loc_cap.dat"])
                    except:
                        print("Can't run OPL model (sp_loc_cap)")
                        exit(1)

                    for safety_margin in safety_margin_range:
                        f = open('SP_cover.dat', "w")
                        f.write(f"n = {n};\n")
                        f.write(f"m = {m};\n")
                        f.write(f"C_num = {C_num};\n")
                        f.write(f"C = {C};\n")
                        f.write(f"h = {h};\n")
                        f.write(f"mu = {mu};\n")
                        f.write(f"p = {p};\n")
                        f.write(f"safety_margin = {safety_margin};\n")
                        f.write(f"Fd = {list(Fd.values())};\n")
                        f.write(f"Df = {list(Df.values())};\n")
                        f.close()

                        f = open(file_name, "a")
                        f.write(f",{safety_margin},")
                        f.close()

                        try:
                            subprocess.run(["oplrun", "sp_cover.mod", "SP_cover.dat"])
                        except:
                            print("Can't run OPL model (sp_cover)")
                            exit(1)

                        " Read solution of the cover model"
                        f = open("out_cover.txt","r")
                        ez = [ float(s) for s in f.readlines()]
                        SP_demand = np.array(ez[:m])
                        SP_type = np.array([int(i) for i in ez[m:]])

                        " estimate number of rjections in solution of the covering model "
                        Rej = 0
                        for i in range(m):
                            if SP_type[i]:
                                Rej += np.interp(SP_demand[i]/safety_margin, Lambda[(SP_type[i]-1),:], R[(SP_type[i]-1),:])

                        f = open(file_name, "a")
                        f.write(f",{Rej}")
                        f.close()

                    f = open(file_name, "a")
                    f.write(f"\n")
                    f.close()
