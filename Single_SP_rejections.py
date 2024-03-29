import numpy as np
from scipy.stats import binom, poisson
import time
import pickle

"""  Ths script accompanied the paper "The service points location and capacity problem" by Tal Raviv (2023)

    It is aimed to evaluate the rejection function.  A similar script "Single_SP_postponement" is used to evaluate
    the postponement function 

    It has two purposes 
    1) Run the experiment reported in Section Appendix
    2) To be used as package with the functions "markov_chain_sp(C, lamb, p)" documented below
        "sim_sp(C, lamb, pd, N=1000000, Warmup=100, rand_seed=0)"
    
    Copyright Tal Raviv (2023). Feel free to use it but please let my know if you do  (talraviv@tau.ac.il)
    Last updated: Feb 17, 2023
"""


"""  markov_chain_sp(C, lamb, p)
        C - Station capacity
        lamb - supply rate
        p - pickup probability 
    
    The function returns 
        R - the mean number of rejections per period  
        CPU_time - the computation time
"""
def markov_chain_sp(C, lamb, p):

    if lamb == 0:
        return 0,0 # don't waste time

    my_inf = int((C+lamb)*2)  # large enough to ignore the tail of the Poisson distribution

    start_time = time.time()
    P = np.zeros((2*C+2, 2*C+2))  # Initialized one step transition probability matrix
    """ Populate the transition probability matrix """
    for i in range(C+1):
        for j in range(i+1):
            P[i, C+1+j] = binom.pmf(i-j, i, p)
            if i < C:
                P[C+1+j,i] = poisson.pmf(i-j,lamb)
            else:
                P[C+1+j,C] = 1- poisson.cdf(C-j-1, lamb)

    ssp = np.linalg.matrix_power(P, 512 )[-1, (C+1):]  # Calculate the before replenishment "steady state" probabilities

    """ Calculate R (not so pythonic) """
    R = 0
    for j in range(C+1):
        for k in range(C-j+1, my_inf):
            R += ssp[j]* poisson.pmf(k, lamb) * (k+j-C)

    return R, time.time()-start_time

"""   Simulation model with the emperical time-to-pickup distribution

      sim_sp(C, lamb, pd)
        C - Station capacity
        lamb - supply rate
        pd - time-to-pickup emperical distribution 
        N - number of periods in the simulation (default 1,000,000)
        Warmup - number of warmup periods (default 100)
        rand_seed - (default 0)
        
    The function returns 
        avg - the estimated mean number of rejections per period  
        hci99, hci95 - half of the width of the 99 and 95 confidence intervals 
        CPU_time - the computation time
"""
def sim_sp(C, lamb, pd, N=1000000, Warmup=100, rand_seed=0):

    np.random.seed(rand_seed)
    rng = np.random.default_rng()
    start_time = time.time()
    H = len(pd)
    departures = np.zeros(N+H)
    R = np.zeros(N)
    I = 0

    for t in range(N):
        NewArrivals = np.random.poisson(lamb)
        toInsert = min(NewArrivals, C-I)
        I += toInsert
        R[t] = NewArrivals - toInsert
        departures[t:t+H] += rng.multinomial(toInsert,pd)
        I -= departures[t]

    """  report sumulation statistics """
    CPU_time =time.time()-start_time
    avg = np.average(R[Warmup:])
    std = np.std(R[Warmup:])
    return avg, std, CPU_time

def erlang_b(E, m):
    inv_b = 1.0
    for j in range(1,m+1):
        inv_b = 1.0 + inv_b * j / E
    return 1.0 / inv_b


"""  Run the small expierment of Section 3.3 in the paper - June 2022 version"""
results = []
SN = 1000000 # number of periods in the simulation (warmup is left at it default, 100)
if __name__ == '__main__':
    Z99 = 2.576  # Z values for confidence interval (assuming N is large)
    Z95 = 1.96
    print("H, p, C, rho, lambda, R_DTMC, cpu_time_DTMC, R_sim, std_sim, hci99, hci95, cpu_time_sim, err_b_approx")
    for pd in [np.array([0.5, 0.2, 0.3]), np.array([0.4, 0.2, 0.1, 0.08, 0.05, 0.02, 0.15])]:# time-to-pickup distribution
        H = len(pd)  # max pickup  periods
        p = 1/np.dot(pd, np.array(range(1,H+1)))   # pickup probability  (geometric time-to-pickup)
        for C in [20, 50, 100]:   # SP capacity
            for rho in [0.9, 0.95, 0.99, 1, 1.01, 1.05, 1.1]:

                lamb = rho*C*p   # supply rate
                err_b_approx = erlang_b(lamb/p, C) * lamb  # not working
                #print(f"Running experiment with H={H}, p={p:.4f}, C={C}, rho={rho}, lambda={lamb:.4f}")
                R_markov, cpu_time_markov =  markov_chain_sp(C, lamb,p)
                avg, std, cpu_time_sim = sim_sp(C, lamb, pd, N=SN)
                hci99 = Z99 * std / np.sqrt(SN-100)
                hci95 = Z95 * std / np.sqrt(SN-100)
                print(f"{H}, {p}, {C}, {rho}, {lamb}, {R_markov}, {cpu_time_markov}, {avg}, {std}, {hci99}, {hci95}, {cpu_time_sim}, {err_b_approx}")
                results.append([H, p, C, rho, lamb, R_markov, cpu_time_markov, avg, std, hci99, hci95, cpu_time_sim, err_b_approx])

pickle.dump(results, open("results_rejection.p", "wb"))
