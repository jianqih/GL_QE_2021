# Algorithm

Step 1: Guess aggregate capital $K0$ and aggregate labor $L0$. 

Step 2: Given $K0$ and $L0$, compute prices $R,w$.

Step 3: Given prices $R,w$ solve HH problems and get policy functions for savings, occupational, and schooling choices. 

Step 4: Given the solution to the HH problem and given the exogenous initial distribution of assets, compute the aggregate saving $K1$ and labor $L1$.

Step 5: $dev V = [K1-K0,L1-L0]$ If devV is above the tolerance level, update guess for $K$ and $L$, and go back to step 2.



HHSimulation_olgm

The basic idea is

(1) Populate a set of households

(2) Households go through a sequence of initial wealth distribution and ability in CovD.

(3) Compute capital holdings of these households based on policy function 



I assume that the occupational path can not transfer with each other freely. i.e. if one person chooses one job after graduation, he cannot change his path even in his ending life. He will continue to work at the 



So in one stage, like retirement period, we only can update the value function via the single career path. i.e. 
$$
vR(states;t)=max \quad u(\cdot)+\beta\cdot vR(states';t+1)
$$






