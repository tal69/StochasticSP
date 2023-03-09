import numpy as np
from scipy.stats import binom, poisson
import time
import MMc
import pickle


"""  Ths script accompanied the paper "The service points location and capacity problem" by Tal Raviv (2023)
    It is used to create the results presented in Appendix B

    It is aimed to evaluate the postponement function.  A similar script "Single_SP" is used to evaluate
    the rejection function 

"""


"""  markov_chain_sp_post(C, lamb, p)
        C - Station capacity
        lamb - supply rate
        p - pickup probability 
    
    The function returns 
        R - the mean number of rejections per period  
        CPU_time - the computation time
        
    NOT WORKING YET    
"""
def markov_chain_sp_post(C, lamb, p, state_factor=5):

    if lamb == 0:
        return 0,0 # don't waste time

    #my_inf = int((C+lamb)*2)  # large enough to ignore the tail of the Poisson distribution
    CC = state_factor*C # the maximal number of parcels in the truncated DTMC (including the queue)

    start_time = time.time()
    P = np.zeros((2*(CC+1), 2*(CC+1)))  # Initialized one-step transition probability matrix
    """ Populate the transition probability matrix """
    for i in range(CC+1):   # number of parcels in the system just after the replenishment
        for j in range(i+1):  # number of parcels in the system just before the replenishment
            P[i, CC+1+j] = binom.pmf(i-j, min(i,C), p)
            if i < CC:
                P[CC+1+j,i] = poisson.pmf(i-j,lamb)
            else:
                P[CC+1+j,CC] = 1- poisson.cdf(CC-j-1, lamb)

    ssp = np.linalg.matrix_power(P, 2048 )[0, :(CC+1)]  # Calculate the before replenishment "steady state" probabilities

    """ Calculate D - the postponement function value (not so pythonic) """
    Lq = np.dot(ssp[C+1:], np.array(range(1,CC-C+1)))

    return Lq, time.time()-start_time

"""   Simulation model with the emperical time-to-pickup distribution

      sim_sp(C, lamb, pd)
        C - Station capacity
        lamb - supply rate
        pd - time-to-pickup emperical distribution 
        N - number of periods in the simulation (default 1,000,000)
        days_in_block - number of days in a block for statistics, first block is omitted as a warmup periods (default 100)
        rand_seed - (default 0)
        
    The function returns 
        avg - the estimated mean number of rejections per period  
        hci99, hci95 - half of the width of the 99 and 95 confidence intervals 
        CPU_time - the computation time
"""
def sim_sp_emprical(C, lamb, pd, N=1000000, rand_seed=0, days_in_block = 100):

    np.random.seed(rand_seed)
    rng = np.random.default_rng()
    start_time = time.time()
    H = len(pd)
    departures = np.zeros(N+H)
    Q_stat = np.zeros(N // days_in_block)
    A_stat = np.zeros(N // days_in_block)
    I = 0  # number of parcels in the SP
    Q = 0 # # number of parcels in the SP

    for t in range(N):
        A = np.random.poisson(lamb)
        Q += A  # add new arrival to queue
        toInsert = min(Q, C-I)
        Q -= toInsert
        I += toInsert
        departures[t:t+H] += rng.multinomial(toInsert,pd)

        I -= departures[t]
        Q_stat[t //days_in_block] += Q
        A_stat[t //days_in_block] += A


    """  report simulation statistics """
    CPU_time =time.time()-start_time
    W = Q_stat / A_stat
    avg = np.average(W[1:])
    std = np.std(W[1:])
    return avg, std, CPU_time


"""   Simulation model with the geometric time-to-pickup distribution

      sim_sp(C, lamb, pd)
        C - Station capacity
        lamb - supply rate
        p - pickup probability at each period 
        N - number of periods in the simulation (default 1,000,000)
        days_in_block - number of days in a block for statistics, first block is omitted as a warmup periods (default 100)
        rand_seed - (default 0)
        
    The function returns 
        avg - the estimated mean number of rejections per period  
        hci99, hci95 - half of the width of the 99 and 95 confidence intervals 
        CPU_time - the computation time
"""
def sim_sp_geometric(C, lamb, p, N=1000000, rand_seed=0, days_in_block = 100):

    np.random.seed(rand_seed)
    rng = np.random.default_rng()
    start_time = time.time()
    H = len(pd)

    Q_stat = np.zeros(N // days_in_block)
    A_stat = np.zeros(N // days_in_block)
    I = 0  # number of parcels in the SP
    Q = 0 # # number of parcels in the SP

    for t in range(N):
        A = np.random.poisson(lamb)
        Q += A  # add new arrival to queue
        toInsert = min(Q, C-I)
        Q -= toInsert
        I += toInsert
        I -= rng.binomial(min(C,I),p)
        Q_stat[t //days_in_block] += Q
        A_stat[t //days_in_block] += A


    """  report simulation statistics """
    CPU_time =time.time()-start_time
    W = Q_stat / A_stat
    avg = np.average(W[1:])
    std = np.std(W[1:])
    return avg, std, CPU_time


"""  Run the experiment of Appendix B in the paper - Feb 2023 version"""

if __name__ == '__main__':


    Z99 = 2.576  # Z values for confidence interval (assuming N is large)
    Z95 = 1.96

    results = []
    print("{H}, {p}, {C}, {rho}, {lamb}, {avg}, {std}, {hci99}, {hci95}, {cpu_time_sim}, {avg_geo}, {hci99_geo}, {hci95_geo}, {cpu_time_sim_geo},{Lq/lamb}, {cpu_time_markov}, {thisQueue.getAvgQueueTime()}")
    for pd in [np.array([0.5, 0.2, 0.3]), np.array([0.4, 0.2, 0.1, 0.08, 0.05, 0.02, 0.15])]:# time-to-pickup distribution
        H = len(pd)  # max pickup  periods
        mean_pd = np.dot(pd, np.array(range(1,H+1)))
        mean_sqr_pd = np.dot(pd, np.square(np.array(range(1,H+1))))
        var_pd = mean_sqr_pd - mean_pd**2
        cv_pd = np.sqrt(var_pd)/mean_pd
        p = 1/mean_pd   # pickup probability  (geometric time-to-pickup)
        for C in [20, 50, 100]:   # SP capacity
            for rho in [0.8,0.9, 0.95, 0.98]:
                lamb = rho*C*p   # supply rate
                thisQueue = MMc.MMcQueue(lamb,p,C)
                W_markov = thisQueue.getAvgQueueTime()
                W_G = 0.5*(cv_pd**2 + 1)*W_markov  # applying Kingman's law of congestion
                #print(f"Running experiment with H={H}, p={p:.4f}, C={C}, rho={rho}, lambda={lamb:.4f}")
                Lq, cpu_time_markov =  markov_chain_sp_post(C, lamb,p)
                avg, std, cpu_time_sim = sim_sp_emprical(C, lamb, pd,N=10000)
                hci99 = Z99 * std / np.sqrt(9999)
                hci95 = Z95 * std / np.sqrt(9999)

                avg_geo, std_geo, cpu_time_sim_geo = sim_sp_geometric(C, lamb, p,N=10000)
                hci99_geo = Z99 * std_geo / np.sqrt(9999)
                hci95_geo = Z95 * std_geo / np.sqrt(9999)


                print(f"{H}, {p}, {C}, {rho}, {lamb}, {avg}, {std}, {hci99}, {hci95}, {cpu_time_sim}, {avg_geo}, {hci99_geo}, {hci95_geo}, {cpu_time_sim_geo},{Lq/lamb}, {cpu_time_markov}, {W_G}")
                results.append([H, p, C, rho, lamb, avg, std, hci95, cpu_time_sim, avg_geo, hci95_geo, cpu_time_sim_geo, {Lq / lamb}, cpu_time_markov, W_G])

pickle.dump(results, open("results_post.p", "wb"))
