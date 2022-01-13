
# pulsating_active_matter
source codes for the manuscript *Pulsating Active Matter*
## Comments on the files
- Param_order.cpp: The C++ code for particle-based simulations, the results of whom were used for measuring phase diagrams.
- simu.sh: The shell script inputting simulation parameters for the C++ code Param_order.cpp
- ep_inc_Param_order.c: The C code for particle-based simulations with **_increasing_** <img src="https://latex.codecogs.com/svg.image?\varepsilon&space;" title="\varepsilon " />
- ep_dec_Param_order.c: The C code for particle-based simulations with **_decreasing_** <img src="https://latex.codecogs.com/svg.image?\varepsilon&space;" title="\varepsilon " />
- ep_inc_simu.sh: The shell script inputting simulation parameters for the C code ep_inc_Param_order.cpp
- ep_dec_simu.sh: The shell script inputting simulation parameters for the C code ep_dec_Param_order.cpp
- movie.py: The Python 3.7 script for trajectory animation
