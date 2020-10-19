##################################################

Paper: 
A mixed hidden Markov model for multivariate monotone disease processes in the presence of measurement errors

Authors:
Lizbeth Naranjo (1), Emmanuel Lesaffre (2), Carlos J. Pérez (3)

Journal:
Statistical Modelling

(1) Departamento de Matemáticas, Facultad de Ciencias, Universidad Nacional Autónoma de México, Ciudad de México, Mexico.
(2) L-BioStat, School of Public Health, KU Leuven, Leuven, Belgium.
(3) Departamento de Matemáticas, Facultad de Veterinaria, Universidad de Extremadura,  Cáceres, Spain.

##################################################

Instructions to run the codes in R and JAGS are provided. 
The codes are applied to obtain a similar analysis as in Section 5 ‘Simulation example’, but with simulated data and without cross-validation. 

##################################################
FILES 

The file ‘MisMHMM.R’ contains the R code. The JAGS code is run from this R file.

The files ‘MisMHMM.jag' and ‘MisMHMM.bug' contains the JAGS model. 

##################################################

To run the files, do the following.
 
1.- Download JAGS from www.mcmc-jags.sourceforge.net/

2.- Install the packages necessary to run the R file. 
These are indicated in the R file. 

3.- Change the address indicated in ‘setwd()’. 
setwd("HERE"). This is the address where the file ‘MisMHMM.bug’ is in.

4.- Run the R file, to simulate the data and to estimate the parameters of the MHMM model.

##################################################
