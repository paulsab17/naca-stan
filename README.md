# NACA-Stan-Modeling

Stan model of a thiol exchange reaction between N-acetylcysteine amide (NACA) and 4,4-aldrithiol (DTP), resulting in an anion of 4-thiopyridine (TP). Additionally, the TP undergoes photobleaching. Absorbance at 325nm is used to calculate [TP] continuously during the reaction, which is performed at several concentrations of the denaturant Guanadinium chloride (Gdn) This analysis uses Monte Carlo sampling and Bayesian statistical packages to estimate the kinetic parameters of the reaction and of photobleaching (more specifically, their values at [Gdn] = 0 and their respective m values).  

nacaAnalysis.Rmd    -   R markdown notebook containing R code for analysis  
testModel.stan      -   stan file describing model to define log probability density function  
reactionDiagram.pdf -   diagram of reaction of interest