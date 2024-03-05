% After running this part an error occurs, but all variables are loaded
run('/Users/sam/Documents/MATLAB/code_2/matlab/main.m');
%%
% patient(404).test = [2867,2502,1731,3568,437,3426,3978,1122,1019,55,3118,987];
% patient(404).test = [2867,2502,1731,3568,437,3426,3978,1122,1019,55,3118,987,3333,3714];
patient(404).test = [2867,2502,1731,3568,437,3426,3978,1122,1019,55,3118,987,3333,3714,2702];
% patient(404).tvalue = ["true","true","true","false","true","true","false","true","true","false","false","false"];
% patient(404).tvalue = ["true","true","true","false","true","true","false","true","true","false","false","false","true","false"];
patient(404).tvalue = ["true","true","true","false","true","true","false","true","true","false","false","false","true","false","true"];
i = 404; 
quickscore_run = 1; 
run("a_quickscore_Alisa.m"); 


