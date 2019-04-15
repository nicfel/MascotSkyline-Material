# Inference of predictor contributions from phylogenetic trees

To recreate the analyses:
1. Run the matlab file _createMasterStepwise.m_, which will create Master input xml's and will output the true values of the coefficients and predictors into files.
2. Run the Master xml files in Beast2 with Master installed.
3. Run the matlab file _createMascotStepwise.m_, which will create Mascot input xml's.
4. Run the Mascot xml file in Beast2 with Mascot installed.
5. Copy the Mascot log files into the out folder.
6. Run _plotScalers.R_ in R to recreate the plots
