# N-body-integrator-comparison

## Overview 
N-body gravitational simulator comparing Euler and leapfrog integration, showing ~1000x lower relative energy drift for leapfrog over 10-200 year solar system simulations. 

## Physics 
This project requires the use of many equations regarding to physics, most notably being Newton's law of gravitation and pairwise force calculations. Numerically speaking, division by zero leads to an undefined result and to ensure this never happens, I've added a little padding to the distance used in calculation. By adding a value so small, it's almost negligible, not affecting the calculations, but ensuring there's no division by zero. The term "integration," in this context, is referring to the process of turning acceleration into velocity into position, over discrete timesteps. 

## Euler vs Leapfrog
Some context. Euler integration is the simplest numerical technique used to approximate solutions to ordinary differential equations. This can lead to some inaccuracy in the final result, as we're not really updating values as much as we should. By utilizing leapfrog integration, we update velocity in kicks (change in velocity due to acceleration), before finding acceleration, instead of using its the previous velocity calculation. Leapfrog integration greatly decreases energy drift, which is the deviation of total system energy, producing a more accurate orbital path.  

## Results
As a result of switching from Euler integration to Leapfrog integration, I have reduced energy drift by a factor of ~1000x over the course of a 10-year, and 200-year, 6-body simulation. Both  integrators show bounded (non-secular) oscillation at these timescales and step sizes. For this experiment, I used a simulation time, dt (how long each step is), paired with the actual steps (how many dt), to accurately simulate the orbital path of a 6 body system. The dt value I used was 86400, which is the amount of seconds in 1 day, along with 3650 steps for 10 years, and 73000 steps for 200 years. 

## How to run
​```
pip install numpy matplotlib 
python gravity_sim_euler.py 
python gravity_sim_leapfrog.py
​```

## Planned improvements
For future improvements, I would like to consolidate both files into one to really highlight the difference in integration technique, and supporting results. I would also like to add energy/angular-momentum conservation diagnostics, as furthermore evidence of Leapfrog's accuracy over Euler. 
