/*********************************************
 * OPL 12.9.0.0 Model
 * Author: talra
 * Creation Date: Jun 22, 2022 at 2:20:59 PM
 *********************************************/
int n = ...;       // Number of demand points
int m = ...;       // Number of candidate SP locations
int C_num = ...;  // Number of possible SP sizes
float alpha = ...;  // Weight of rejections in the objective function

range DemandPoints = 0..n-1;  // set of demand points
range SP_loc = 0..m-1;    // set of candidate locations
range SP_cap = 1..C_num;  //  set of possible SP Capacities

float mu[DemandPoints] = ...;
float h[SP_loc, SP_cap] = ...;
int Lambda_points = ...;
float Lambda[SP_cap,1..Lambda_points] = ...;
float R[SP_cap, 1..Lambda_points] = ...;
{int} Fd[DemandPoints] = ...;
{int} Df[SP_loc] = ...;
float temp;
// Information only 

//float p = ...;
//float r = ...;
//int number_of_xs = ...;

dvar boolean y[SP_loc, SP_cap];
dvar float+ x[SP_loc, DemandPoints];
dvar float+ z[SP_loc, 1..Lambda_points, 1..C_num];

dexpr float SP_cost = sum(s in SP_cap, f in SP_loc)  h[f,s] * y [f,s];
dexpr float rejections = sum(s in SP_cap, f in SP_loc, k in 1..Lambda_points) z[f,k,s] * R[s,k];
dexpr int Opened_SPs[s in SP_cap] = sum(f in SP_loc) y[f,s]; 

execute {cplex.tilim = 600;
	var before = new Date();
	temp = before.getTime();
}

minimize SP_cost + alpha*rejections;
subject to {
	forall ( f in SP_loc) sum( s in 1..C_num) y[f,s] <= 1;  // Constraint (8)
  	forall(d in DemandPoints) sum(f in Fd[d]) x[f,d] == mu[d];   // Constraint (9)
  	forall(f in SP_loc) sum(k in 1..Lambda_points, s in 1..C_num) z[f,k,s] * Lambda[s,k] 
  		== sum(d in Df[f]) x[f,d];   // Constraint (10)
  		
  	forall(f in SP_loc, s in SP_cap) sum(k in 1..Lambda_points) z[f,k,s] == y[f,s];  // Constraint (11)
}  


execute 
{	
	var after = new Date();
	var f = new IloOplOutputFile("Results_SP.csv", true); // result file
	f.write(cplex.getCplexStatus(),","); // exit status 1-ok; 11 - stop due to time limit
	f.write(cplex.getObjValue(),",");  // obj_function
	f.write(cplex.getBestObjValue(),","); // LB
	f.write(after.getTime()-temp,",")
	f.write(SP_cost);  // obj func components
	for (var s in SP_cap) f.write(",",Opened_SPs[s]);
	f.write(",",rejections);  // obj func components
	// Line feed will be added by the Pyhton script
	f.close();
}	