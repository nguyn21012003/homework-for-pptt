# sbesolidphysics
#DECLARATION <br />
This code is done by the main contribution from the theoretical group B38 (at now we will call by TTCQV + D group). <br />
The copy of this code with the other version is NO ACCCEPTED. <br />

#AIM <br />
The programing solves sbe problem (system of ode) and generates the target imange via Rugga Kutta 4th order method.  <br />

#FUNCTIONS <br />
Using the following command to run this code <br />

With make <br />
 - Ensuring that you have already installed "make" <br />
 - "make run" to run the full code. <br />

Without make <br />
 You can run it yourself by the following command
 - Complide and run file sbe.c or sbe_v1.c by "gcc sbe.c -o sbe.out -lm" (or sbe_v1.c). <br />
 - Run "./sbe.out". <br />
 - After checked the file ".txt" exist, you can run the "python3 plot.py" to generate the two images. <br />
 - Done. <br />

#RESULTS <br />
The output file are two images (the distribution and polarization function in term of time) and one .txt data file <br />
 - Plotting style can be changed in plot.py. <br />
Please put full code in one folder to make code clear run. <br />

