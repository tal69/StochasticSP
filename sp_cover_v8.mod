/*********************************************
 * OPL  Model
 * Author: talra
 * Creation Date: Feb 21, 2023 
 *********************************************/
int time_limit  = ...;
int n = ...;       // Number of demand points
int m = ...;       // Number of candidate SP locations
int C_num = ...;  // Number of possible SP sizes


float alpha = ...;  // Weight of rejections in the objective function
//float r = ...;  // service radius
float p = ...;   // mean time in SP
float safety_margin = ...;
range DPs = 0..n-1;  // set of demand points
range SPs = 0..m-1;    // set of candidate locations



range SP_cap = 1..C_num;  //  set of possible SP Capacities
int C[SP_cap] = ...;

float mu[DPs] = ...;
float h[SPs, SP_cap] = ...;    // set up costs
int Lambda_points = ...;           // number of discretization points   - only for compatability here
float Lambda[SP_cap,1..Lambda_points] = ...;// discretization of the rejection function (x-axis)  - only for compatability here
float R[SP_cap, 1..Lambda_points] = ...; // discretization of the rejection function (y-axis)  - only for compatability here

tuple T_DP2SP {
  key int dp;
  key int sp;
}

{T_DP2SP} DP2SP = ...;

{int}farther_SPs[DP2SP] = ...;



float temp;

dvar boolean y[SPs, SP_cap];
dvar float+ x[ DP2SP];  
dvar float+ q[SPs];

dexpr float SP_cost = sum(s in SP_cap, f in SPs)  h[f,s] * y [f,s];
dexpr float rejections = sum(f in SPs) q[f];
dexpr int Opened_SPs[s in SP_cap] = sum(f in SPs) y[f,s];
dexpr int SP_sol[f in SPs] = sum(s in SP_cap) s* y[f,s]; 
dexpr float demand_SP[f in SPs] = sum(i in DP2SP : i.sp == f) x[i];

execute {cplex.tilim = time_limit;
	var before = new Date();
	temp = before.getTime();
}

minimize SP_cost + alpha * rejections;
subject to {


	forall(f in SPs) sum(i in DP2SP : i.sp == f) x[i] <= sum(s in SP_cap) p*C[s]* y[f,s] + q[f]; 
	// enough capacity		

	forall(d in DPs) sum(i in DP2SP : i.dp == d)  sum( s in SP_cap)  y[i.sp,s] >= 1; // coverage constraint
	
	forall(d in DPs) sum(i in DP2SP: i.dp == d) x[i] == safety_margin*mu[d];   // Demand satisfaction
	
	forall ( f in SPs) sum( s in SP_cap) y[f,s] <= 1;  // At most one APL at each SP
  	
  	// supply only from cloesest SP
	forall(i in DP2SP )  
		sum(f in farther_SPs[i]) x[<i.dp,f>] <= mu[i.dp]*(1- sum(s in SP_cap) y[i.sp,s]);
	
}  



execute 
{
	
	var after = new Date();
	
	

	var f = new IloOplOutputFile("Results_SP_v8.csv", true); // result file
	f.write(",",safety_margin,",",cplex.getCplexStatus(),","); // exit status 1-ok; 11 - stop due to time limit
	f.write(cplex.getObjValue(),",");  // obj_function
	f.write(cplex.getBestObjValue(),","); // LB
	f.write(after.getTime()-temp,",")
	f.write(SP_cost,",",rejections);  // obj func components
	for (var s in SP_cap) f.write(",",Opened_SPs[s]);
	// Line feed will be added by the Pyhton script
	f.close(); 

	f = new IloOplOutputFile("out_cover.txt"); // result file
	for(var i in SPs) f.writeln(SP_sol[i]);
	for(var i in SPs) f.writeln(demand_SP[i]);
	f.close();	
	
	
}	
