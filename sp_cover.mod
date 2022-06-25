/*********************************************
 * OPL 12.9.0.0 Model
 * Author: talra
 * Creation Date: Jun 23, 2022 at 2:20:59 PM
 *********************************************/
int n = ...;       // Number of demand points
int m = ...;       // Number of candidate SP locations
int C_num = ...;  // Number of possible SP sizes

range DemandPoints = 0..n-1;  // set of demand points
range SP_loc = 0..m-1;    // set of candidate locations
range SP_cap = 1..C_num;  //  set of possible SP Capacities
int C[SP_cap] =...;
float p = ...;
float mu[DemandPoints] = ...;
float h[SP_loc, SP_cap] = ...;
{int} Fd[DemandPoints] = ...;                
{int} Df[SP_loc] = ...;
float safety_margin = ...;
float temp;

dvar boolean y[SP_loc, SP_cap];
dvar float+ x[SP_loc, DemandPoints];

dexpr float SP_cost = sum(s in SP_cap, f in SP_loc)  h[f,s] * y [f,s];
dexpr int Opened_SPs[s in SP_cap] = sum(f in SP_loc) y[f,s]; 
dexpr float demand_SP[f in SP_loc] = sum(d in Df[f]) x[f,d];
dexpr int SP_sol[f in SP_loc] = sum(s in SP_cap) s* y[f,s]; 
execute {cplex.tilim = 600;
	var before = new Date();
	temp = before.getTime();
}

minimize SP_cost;
subject to {
	forall ( f in SP_loc) sum( s in SP_cap) y[f,s] <= 1;  // At most one APL at each SP
  	forall(d in DemandPoints) sum(f in Fd[d]) x[f,d] == safety_margin*mu[d];   // Demand satisfaction
  	forall(f in SP_loc) sum(d in Df[f]) x[f,d] <= sum(s in SP_cap) p*C[s]* y[f,s];   // Constraint (16)		
}  


execute 
{
	
	var after = new Date();

	var f = new IloOplOutputFile("Results_SP.csv", true); // result file
	f.write(cplex.getCplexStatus(),","); // exit status 1-ok; 11 - stop due to time limit
	f.write(cplex.getObjValue(),",");  // obj_function
	f.write(cplex.getBestObjValue(),","); // LB
	f.write(after.getTime()-temp)
	for (var s in SP_cap) f.write(",",Opened_SPs[s]);
	// Line feed will be added by the Pyhton script
	f.close();
	
	f = new IloOplOutputFile("out_cover.txt"); // result file
	for(var i in SP_loc) f.writeln(demand_SP[i]);
	for(var i in SP_loc) f.writeln(SP_sol[i]);
	f.close();
	
}	