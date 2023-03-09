/*********************************************
 * OPL 12.9.0.0 Model
 * Author: talra
 * Creation Date: Jun 22, 2022 at 2:20:59 PM
 *********************************************/
int time_limit  = ...;
int n = ...;       // Number of demand points
int m = ...;       // Number of candidate SP locations
int C_num = ...;  // Number of possible SP sizes
float alpha = ...;  // Weight of rejections in the objective function
float r = ...;  // service radius

range DemandPoints = 0..n-1;  // set of demand points
range SP_loc = 0..m-1;    // set of candidate locations
range SP_cap = 1..C_num;  //  set of possible SP Capacities

float mu[DemandPoints] = ...;
float h[SP_loc, SP_cap] = ...;
int Lambda_points = ...;
float Lambda[SP_cap,1..Lambda_points] = ...;
float R[SP_cap, 1..Lambda_points] = ...;
float SP_xy[SP_loc, 1..2] = ...;
float DP_xy[DemandPoints, 1..2] = ...;

float temp;

float dist[d in DemandPoints, f in SP_loc] = sqrt(sum(k in 1..2) pow(SP_xy[f,k]- DP_xy[d,k],2));

dvar boolean y[SP_loc, SP_cap];
dvar float+ x[ DemandPoints, SP_loc];   // todo - define only on the relevant tuples where dist(d,f)< r
dvar float+ z[SP_loc, 1..Lambda_points, 1..C_num];

dexpr float SP_cost = sum(s in SP_cap, f in SP_loc)  h[f,s] * y [f,s];
dexpr float rejections = sum(s in SP_cap, f in SP_loc, k in 1..Lambda_points) z[f,k,s] * R[s,k];
dexpr int Opened_SPs[s in SP_cap] = sum(f in SP_loc) y[f,s]; 
//dexpr int SP_sol[f in SP_loc] = sum(s in SP_cap) s* y[f,s]; 


execute {cplex.tilim = time_limit;
	var before = new Date();
	temp = before.getTime();
}

minimize SP_cost + alpha*rejections;
subject to {
	forall ( f in SP_loc) sum( s in 1..C_num) y[f,s] <= 1;  // At most one capacity for each candiate location
  	forall(d in DemandPoints) sum(f in SP_loc: dist[d,f]<r) x[d,f] == mu[d];   // Demand satisfaction
  	forall(f in SP_loc) sum(k in 1..Lambda_points, s in 1..C_num) z[f,k,s] * Lambda[s,k] 
  		== sum(d in DemandPoints : dist[d,f] < r) x[d,f];   // Assign correct value to Z
  		
  	forall(f in SP_loc, s in SP_cap) sum(k in 1..Lambda_points) z[f,k,s] == y[f,s];  // Sum of relevant z's equal 1, others 0

	forall(d in DemandPoints, f in SP_loc : dist[d,f] < r )  forall(ff in SP_loc: dist[d,ff] < dist[d,f]) 
	 	 sum(fff in SP_loc : dist[d,fff] >= dist[d,f] ) x[d,fff] <=
		 mu[d] * (1- sum(s in SP_cap) y[ff,s]);

}  


execute 
{	
	var after = new Date();
	var f = new IloOplOutputFile("Results_SP_v5.csv", true); // result file
	f.write(",",cplex.getCplexStatus(),","); // exit status 1-ok; 11 - stop due to time limit
	f.write(cplex.getObjValue(),",");  // obj_function
	f.write(cplex.getBestObjValue(),","); // LB
	f.write(after.getTime()-temp,",")
	f.write(SP_cost,",",rejections);  // obj func components
	for (var s in SP_cap) f.write(",",Opened_SPs[s]);
	// Line feed will be added by the Pyhton script
	f.close();
	
}	
