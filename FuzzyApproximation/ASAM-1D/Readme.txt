Author (version 2.1) : Olaoluwa Adigun
Author (version 2.0) : Osonde Osoba
This is an update of version 2.0 with additional features.

Version (2.0): Trains 6 SAMs to approximate a user-specified 1-dimensional function using samples of the objective function.
Version (2.1): The added features are:
              --- Computing the conditonal variance for the Fuzzy approximation.
              --- Fuzzy aproximation for noisy functions 
         I
How to run:
===========
1) - The program takes a function f(x) and finds its fuzzy approximation. You can specify the desired function.
x	f(x)
for example:
5.56565657	0.0929060987	
I have included the sample "FxnGen.dat" file for you to compare.

2) - The program outputs 6 folders: Gauss, Sinc, Cauchy, etc. These contain the ASAM approximation details for the three different set functions. You should tune the number of rules, number of iterations, learning rates, and initializations to get better performance on your approximations. It also outputs the "Errors.dat" file. This is a log of MSEs for all the ASAMs at different points in the adaptation. "OutFxn.dat" is just a sanity-check procedure. 

3) There is a make file for running the code. Just type "make ASAM1D" on the CMD terminal.


Caveats:
========
1) - The program erases outputs from previous runs every time you run it. You might need cygwin on Windows to deal for folder manipulation. 
