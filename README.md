# NumTechFluidDyn

# PROGRESS:
- Session 1-4: finished
- Session 5 (Method of manufactured solutions): bugs in implementation. 
	I implemented the method of manufactured solutions. We now have 3 fields in ExampleCaseSession5: 
		- U_sol is the field with the exact analytic results
		- U is the solution field computed by the code
		- U_diff is the error field
	The source term computed in the slides contains errors, Ishaan forgot to include the viscosity v in one of the terms
	I implemented this already correctly in the code.
	When including the contribution in the bdata matrix, don't forget that the source term gets integrated over the cell volume. 
	This means you should multiply it with cVol.
	RESULTS: run ExampleCaseSession5.m to see the results (with solver exampleSolver5). 
	It seems like the boundary condition is implemented correctly, but there is something wrong (look at the plots, it will be clear).
	However, I think the error is not in the implementation of the MMS, but in the body of our code. I assume it will have something to do with the coupling of the u and v component of the velocity 
	(this would explain why we got the exact results in the previous sessions: the v component was 0 there).