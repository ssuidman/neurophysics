% After running this part an error occurs, but all variables are loaded
run('/Users/sam/Documents/MATLAB/code_2/matlab/main.m');
%%
% Then this part can be ran in which a patient 999 is created, the
% quickscore algorithm is run and the variables are saved in a mat-file 

% Find cases: 
% patient(406).tvalue{2} = 'false' % this case has a mistake in it
% for i=1:1000; k=0; for j=1:size(patient(i).tvalue,2); if ~contains(['normal','false'],patient(i).tvalue{j});k=k+1; end; end; if k~=0; fprintf("i=%i m=%i\n",i,k); end; end;


% patient(998).test = [55,437,987,1019,1122,1420,1731,2867,3118,3235,3284,3426,3568,3812,3978,2687,2689,2690,2691,2692,2699,2702,2704,2706];
% patient(998).tvalue = ["true","true","true","true","true","true","true","true","true","true","true","true","false","false","true",   "true","true","true","true","true","true","true","true","true"];
% patient(998).age = 40; 
% patient(998).gender = 'F';
% i = 998; 

% patient(997).test = [55,437,987,1019,1122,1420,1731,2867,3118,3235,3284,3426,3568,3812,3978,2687,2689];
% patient(997).tvalue = ["true","true","true","true","true","true","true","true","true","true","true","true","false","false","true","true","true"];
% patient(997).age = 40; 
% patient(997).gender = 'F';
% i = 997; 

% Noisy Max case --> patientprev, sens_normal, sens_extreme 
% i = 320; 
% Noisy Max case --> patientprev, sens_normal, sens_extreme 
i = 105; 


% quickscore_run = 0; 
% run('/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/a_quickscore_Alisa.m'); 
% previn = patientprev(:,1)';
% pfmin = sens_extreme';
% pfminneg = sens_normal';
% For the m=22 case of patient 999:
% [pfplus, P_joint, posterior, dt, dt_array] = b_quickscore_Wim(previn,pfmin,pfminneg); 
% dt_array = [0.016340, 0.020360, 0.025232, 0.033703, 0.054161, 0.086878, 0.144844, 0.245685, 0.443075, 0.828686, 1.645363, 3.184414, 6.230765, 12.324080, 24.584273, 49.069105, 97.992280, 197.838568, 395.341900, 790.240980, 1594.027182, 3227.494602]; 
% save("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_997.mat","previn","pfmin","pfminneg","pfplus","P_joint","posterior","dt","dt_array");

% save("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_320.mat","patientprev","sens_normal","sens_medium","sens_extreme");
% save("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_105.mat","patientprev","sens_normal","sens_medium","sens_extreme");
