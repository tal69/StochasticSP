/*********************************************
 * OPL Model
 * Author: Tal Raviv
 * Creation Date: Feb 21, 2023 
 *********************************************/
 
 
int time_limit  = ...;
int n = ...;       // Number of demand points
int m = ...;       // Number of candidate SP locations
int C_num = ...;  // Number of possible SP sizes
float alpha = ...;  // Weight of rejections in the objective function
//float r = ...;  // service radius

range DPs = 0..n-1;  // set of demand points
range SPs = 0..m-1;    // set of candidate locations



range SP_cap = 1..C_num;  //  set of possible SP Capacities



float mu[DPs] = ...;
float h[SPs, SP_cap] = ...;    // set up costs
int Lambda_points = ...;           // number of discretization points
float Lambda[SP_cap,1..Lambda_points] = ...;// discretization of the rejection function (x-axis)
float R[SP_cap, 1..Lambda_points] = ...; // discretization of the rejection function (y-axis)

tuple T_DP2SP {
  key int dp;
  key int sp;
}

{T_DP2SP} DP2SP = ...;

{int}farther_SPs[DP2SP] = ...;



float temp;

dvar boolean y[SPs, SP_cap];
dvar float+ x[ DP2SP];  
dvar float+ z[SPs, 1..Lambda_points, 1..C_num];

dexpr float SP_cost = sum(s in SP_cap, f in SPs)  h[f,s] * y [f,s];
dexpr float rejections = sum(s in SP_cap, f in SPs, k in 1..Lambda_points) z[f,k,s] * R[s,k];
dexpr int Opened_SPs[s in SP_cap] = sum(f in SPs) y[f,s];
dexpr int SP_sol[f in SPs] = sum(s in SP_cap) s* y[f,s]; 
dexpr float demand_SP[f in SPs] = sum(i in DP2SP : i.sp == f) x[i];

execute {cplex.tilim = time_limit;
	var before = new Date();
	temp = before.getTime();
}

minimize SP_cost + alpha*rejections;
subject to {
	forall ( f in SPs) sum( s in 1..C_num) y[f,s] <= 1;  // At most one capacity for each candiate location
  
  	forall(d in DPs) sum(i in DP2SP : i.dp == d) x[i] == mu[d];   // Demand satisfaction
  
  	forall(f in SPs) sum(k in 1..Lambda_points, s in 1..C_num) z[f,k,s] * Lambda[s,k] 
  		== demand_SP[f];   // Assign correct value to Z
  		
  	forall(f in SPs, s in SP_cap) sum(k in 1..Lambda_points) z[f,k,s] == y[f,s];  // Sum of relevant z's equal 1, others 0


	forall(i in DP2SP )  
		sum(f in farther_SPs[i]) x[<i.dp,f>] <= mu[i.dp]*(1- sum(s in SP_cap) y[i.sp,s]);
	 
}  


execute 
{	
	var after = new Date();
	var f = new IloOplOutputFile("Results_SP_v8.csv", true); // result file
	f.write(",",cplex.getCplexStatus(),","); // exit status 1-ok; 11 - stop due to time limit
	f.write(cplex.getObjValue(),",");  // obj_function
	f.write(cplex.getBestObjValue(),","); // LB
	f.write(after.getTime()-temp,",")
	f.write(SP_cost,",",rejections);  // obj func components
	for (var s in SP_cap) f.write(",",Opened_SPs[s]);
	// Line feed will be added by the Pyhton script
	f.close();
	
	f = new IloOplOutputFile("out_our.txt"); // result file
	for(var i in SPs) f.writeln(SP_sol[i]);
	for(var i in SPs) f.writeln(demand_SP[i]);
	f.close();
	
	
}	
