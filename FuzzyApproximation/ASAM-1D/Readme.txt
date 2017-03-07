Author (version 1.1) : Olaoluwa Adigun
Author (version 1.0) : Osonde Osoba
This is an update of version 1.1 with additional features.

Version (1.0): Trains 6 SAMs to approximate a user-specified 1-dimensional function using samples of the objective function.
Version (1.1): The added features are:
              --- Computing the conditonal variance for the Fuzzy aproximation.
              --- Fuzzy aproximation for noisy functions 
         I
How to run:
===========
1) - The program takes function samples from a file "FxnGen.dat". You supply this file in the same directory as the program. The samples are multiple lines in the form:
x	f(x)
for example:
5.56565657	0.0929060987	
I have included the sample "FxnGen.dat" file for you to compare.

2) - The program outputs 6 folders: Gauss, Sinc, Cauchy, etc. These contain the ASAM approximation details for the three different set functions. You should tune the number of rules, number of iterations, learning rates, and initializations to get better performance on your approximations. It also outputs the "Errors.dat" file. This is a log of MSEs for all the ASAMs at different points in the adaptation. "OutFxn.dat" is just a sanity-check procedure. 



Caveats:
========
1) - The program erases outputs from previous runs every time you run it. You might need cygwin on Windows to deal for folder manipulation. It is currently set to use 39 rules and stop after 60000 iterations.
