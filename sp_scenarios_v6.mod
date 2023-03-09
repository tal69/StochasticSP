/*********************************************
 * OPL 12.9.0.0 Model
 * Author: Tal Raviv
 * Creation Date: Jun 27, 2022 
 *********************************************/
int time_limit  = ...;
int n = ...;       // Number of demand points
int m = ...;       // Number of candidate SP locations
int C_num = ...;  // Number of possible SP sizes
int scenario_num = ...; // Number of scenarios
float alpha = ...;  // Weight of rejections in the objective function
float r = ...;  // service radius

range DemandPoints = 0..n-1;  // set of demand points
range SP_loc = 0..m-1;    // set of candidate locations
range SP_cap = 1..C_num;  //  set of possible SP Capacities
range Scenarios = 1..scenario_num; // set of scenarios

int C[SP_cap] =...;
float p = ...;
int mu[DemandPoints,Scenarios] = ...;
float h[SP_loc, SP_cap] = ...;
float SP_xy[SP_loc, 1..2] = ...;
float DP_xy[DemandPoints, 1..2] = ...;
float temp;
float dist[d in DemandPoints, f in SP_loc] = sqrt(sum(k in 1..2) pow(SP_xy[f,k]- DP_xy[d,k],2));

dvar boolean y[SP_loc, SP_cap];
dvar float+ x[ DemandPoints,SP_loc, Scenarios];
dvar float+ q[SP_loc, Scenarios];

dexpr float SP_cost = sum(s in SP_cap, f in SP_loc)  h[f,s] * y [f,s];
dexpr int Opened_SPs[s in SP_cap] = sum(f in SP_loc) y[f,s]; 
dexpr int SP_sol[f in SP_loc] = sum(s in SP_cap) s* y[f,s];
dexpr float rejections = (1/scenario_num) * sum ( r in Scenarios , f in SP_loc) q[f,r];


execute {cplex.tilim = time_limit;
	var before = new Date();
	temp = before.getTime();
}

minimize SP_cost + alpha * rejections;
subject to {

	forall(f in SP_loc, sce in Scenarios) sum(d in DemandPoints : dist[d,f] < r) x[d,f,sce] <= sum(s in SP_cap) p*C[s]* y[f,s]  + q[f,sce];   // Capacity constraint (16)		
	
	forall ( f in SP_loc) sum( s in SP_cap) y[f,s] <= 1;  // At most one APL at each SP
  	
	forall(d in DemandPoints, sce in Scenarios) sum(f in SP_loc : dist[d,f]<r) x[d,f,sce] == mu[d,sce] ;   // Demand satisfaction
	
	forall(d in DemandPoints, f in SP_loc : dist[d,f] < r )  
	 	 sum(ff in SP_loc : dist[d,ff] > dist[d,f] ) sum(r in Scenarios) x[d,ff,r] <=
		 sum( sce in Scenarios) mu[d,sce] * (1- sum(s in SP_cap) y[f,s]);
	
	forall(d in DemandPoints) sum(f in SP_loc : dist[d,f] < r ) sum( s in SP_cap)  y[f,s] >= 1; // coverage constraint
}

execute 
{

	f = new IloOplOutputFile("out_scenarios.txt"); // result file
	for(var i in SP_loc) f.writeln(SP_sol[i]);
	f.close();

	var after = new Date();
	var f = new IloOplOutputFile("Results_SP_v6.csv", true); // result file
	f.write(",", scenario_num)
	f.write(",", cplex.getCplexStatus(),","); // exit status 1-ok; 11 - stop due to time limit
	f.write(cplex.getObjValue(),",");  // obj_function
	f.write(cplex.getBestObjValue(),","); // LB
	f.write(after.getTime()-temp,",")
	f.write(SP_cost,",",rejections);  // obj func components
	for (var s in SP_cap) f.write(",",Opened_SPs[s]);
	f.close();
	
}	
