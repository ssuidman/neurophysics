% After running this part an error occurs, but all variables are loaded
run('/Users/sam/Documents/MATLAB/code_2/matlab/main.m');
%%
for x=1:10
    display(x);
    tests = zeros(1,20); 
    tvalues = strings(1,20); 
    k = 1; 
    while sum(tests ~= 0)~=20
        i = randi([1,4011]); 
        if ismember(5,test(i).type) && ~ismember(i,tests) 
            tests(k) = i;
            if sum(tests ~= 0)<3
                tvalues(k) = "false";
            else
                tvalues(k) = "true";
            end
            k = k+1; 
        end
    end
    [tests,idx] = sort(tests); 
    tvalues = tvalues(idx);
    
    patient(999).test = tests; 
    patient(999).tvalue = tvalues; 
    patient(999).age = 30; 
    patient(999).gender = 'F';
    display(patient(999)); 

    quickscore_run = 0; 
    i = 999; 
    run('/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/a_quickscore_Alisa.m'); 
    previn = patientprev(:,1)';
    pfmin = sens_extreme';
    pfminneg = sens_normal';
    filename = sprintf('/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_%d.mat', x);
%     save(filename, 'previn', 'pfmin', 'pfminneg');
end 

