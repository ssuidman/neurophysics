% After running this part an error occurs, but all variables are loaded
run('/Users/sam/Documents/MATLAB/code_2/matlab/main.m');
%%
% Then this part can be ran in which a patient 999 is created, the
% quickscore algorithm is run and the variables are saved in a mat-file 
patient(999).test = [55,437,987,1019,1122,1420,1731,2867,3118,3235,3284,3426,3568,3812,3978];
patient(999).tvalue = ["true","true","true","true","true","true","true","true","true","true","true","true","false","false","true"];
patient(999).age = 40; 
patient(999).gender = 'F';
i = 999; 
quickscore_run = 1; 
run('/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/a_quickscore_Alisa.m'); 
save("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_999.mat","patientprev","sens_normal","sens_medium","sens_extreme","diagset","pdiag");
