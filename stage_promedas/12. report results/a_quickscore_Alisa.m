fid=fopen('results','w');

% i = 404;
% i = 999;

tic
if (~isempty(patient(i).age))&(length(patient(i).test)==length(patient(i).tvalue))

fprintf(fid,'# Patient N = %d\n',i);
fprintf(fid,'%s %d\n','Age: ',patient(i).age);
switch patient(i).gender,
    case 'M'
        fprintf(fid,'%s %s\n','Gender: ','male');
    case 'F'
        fprintf(fid,'%s %s\n','Gender: ','female');
    end %switch

tests_1=unique(patient(i).test);
if length(tests_1)<length(patient(i).test)
    fprintf(fid,'%s\n','');
    fprintf('ERROR. Test repetition\n');
    fprintf(fid,'ERROR. Test repetition\n');
end

%1 deviding findings into 8 different groups (according to noisy-MAX model)
% creating diagtests - list of diagnoses which we will consider (connected to positive findings)
Below_norm=[]; Below_b=[]; Below_eb=[];
Above_norm=[]; Above_ea=[]; Above_a=[]; 
Binary_true=[]; Binary_false=[]; diagset=[];

diagset=union(diagset,patient(i).diag);

for k=1:length(tests_1)
    t=tests_1(k);
    index=find(patient(i).test==t,1);
    value=patient(i).tvalue{index};
     if ismember (3,test(t).type) %checking for below 
        if strmatch('below',value)
           Below_b(end+1)=t; 
           pparents=find(full(belowSens(:,t))~=0);
           diagset=union(diagset,pparents);
        elseif strmatch('ebelow',value)
           Below_eb(end+1)=t;
%                pparents=find(full(ebelowSens(:,t))~=0);
%                diagset=union(diagset,pparents);
           for j=1:length(parents{t})
               d=parents{t}(j);
               if diagtest(d,t).sens(5)>0
                    diagset=union(diagset,d);
               end
           end
        elseif (strcmp('above',value)||strcmp('eabove',value))||strcmp('normal',value)
           Below_norm(end+1)=t; 
        else
           fprintf('%d  %s  %s\n',t,test(t).name,'ERROR:Inlavid value'); 
        end
    end
    if ismember (2,test(t).type) %checking for above
        if strmatch('above',value)
           Above_a(end+1)=t;
           pparents=find(full(aboveSens(:,t))~=0);
           diagset=union(diagset,pparents);
        elseif strmatch('eabove',value)
           Above_ea(end+1)=t;
%                pparents=find(full(eaboveSens(:,t))~=0);
%                diagset=union(diagset,pparents);
           for j=1:length(parents{t})
               d=parents{t}(j);
               if diagtest(d,t).sens(1)>0
                    diagset=union(diagset,d);
               end
           end
        elseif (strcmp('below',value)||strcmp('ebelow',value))||strcmp('normal',value) 
           Above_norm(end+1)=t;
        else
           fprintf('%d  %s  %s\n',t,test(t).name,'ERROR:Inlavid value'); 
        end
    end
    if test(tests_1(k)).type==5
        if strmatch('true',value)
           Binary_true(end+1)=t;
           pparents=find(full(restSens(:,t))~=0);
           diagset=union(diagset,pparents);
        elseif strmatch('false',value)
           Binary_false(end+1)=t;
        else
           fprintf('%-2d  %-70s %s\n',t,test(t).name,'ERROR:Inlavid value');
        end
    end
end %for k

diagset2 = diagset; 

%2. computing prevalences of considered diagnoses based on patient's age, gender and other conditions: 
patientprev=zeros(length(diagset),1);
for j=1:length(diagset)
    d=diagset(j);       
    % compute patient specific prev for each diag
    switch patient(i).gender,
    case 'M'
        patientprev(j)=prev(d)*malemult(d);
    case 'F'
        patientprev(j)=prev(d)*femalemult(d);
    end %switch

    if ~isempty(agemult(d).low)
        % diag is age dependent
        low_tmp=agemult(d).low;
        high_tmp=agemult(d).high;
        mult_tmp=agemult(d).mult;
        for k=1:length(low_tmp)
            if and(low_tmp(k)<=patient(i).age,patient(i).age<high_tmp(k))
                patientprev(j)=patientprev(j)*mult_tmp(k);
            end
        end %for
    end %if

    % prev multiplication due to externals	
    for k=1:length(patient(i).ext)
        t=patient(i).ext(k);
        if link_ext(d,t)&strmatch('true',patient(i).evalue(k))
            patientprev(j)=patientprev(j)*diagtest(d,t).mult;
        end
    end %for k
end %for j 

%checking that patientprev is smaller than 1: --> because multipliers could
 % make it bigger than 1 
patientprev(patientprev>1)=0.5;
index=find(patientprev>1);
if nnz(index)
    fprintf(fid,'%s\n','');
    fprintf('ERROR. Prior probability >1:  %d\n',diagset(index));
    fprintf(fid,'ERROR. Prior probability >1: diagnosis %d - prob:%f\n',diagset(index),patientprev(index));
end

%setting the prior probability P=1 to known patient's diagnoses:
patientprev=[patientprev,1-patientprev];
fprintf('diagnosis \t prevalence \t 1-prevalence \t sum');fprintf('\n');
for j=1:length(patient(i).diag)
    index=find(diagset==patient(i).diag(j));
    if startsWith(cell2mat(patient(i).dvalue(j)),'tru')
        fprintf('%.10f \t ',[index,patientprev(index,1),patientprev(index,2),sum(patientprev(index,1:2))]);fprintf('\n');
        patientprev(index,2)=0;
        fprintf('%.10f \t ',[index,patientprev(index,1),patientprev(index,2),sum(patientprev(index,1:2))]);fprintf('\n');
    elseif startsWith(cell2mat(patient(i).dvalue(j)),'fals')
        patientprev(index,1)=0;
    else 
        fprintf('%d  %s  %s\n',patient(i).diag(j),cell2mat(diagn(diagset(index))),'ERROR:Inlavid diagnosis value');
    end
end

patientprev3 = patientprev;
index=find(patientprev(:,1)==0);
diagset(index)=[]; 
patientprev(index,:)=[];

%3 creating matrix (pat_sensdt) P(tj=1|di=1) for existing positive and negative findings and
% parent diagnoses:

Normal=[Binary_false, Below_norm, Above_norm]; %negative findings
Medium=[Below_b, Above_a]; %cont tests "above" and "below"
Extreme=[Binary_true, Below_eb, Above_ea]; %%cont tests "eabove", "ebelow" and binary tests "true"

sens_normal=ones(length(diagset),1);
sens_medium=zeros(length(diagset),length(Medium),2);
sens_extreme=zeros(length(diagset),length(Extreme));

for j=1:length(diagset)
    d=diagset(j); 

    %negitive findings( binary and continuous):
    for k=1:length(Binary_false) 
        t=Binary_false(k);
        sens_normal(j)=sens_normal(j)*(1-restSens(d,t)); 
    end
    for k=1: length(Below_norm)
        t=Below_norm(k);
        %sens_normal(j)=sens_normal(j)*(1-belowSens(d,t)); 
        if ~isempty(diagtest(d,t).sens)
            sens_normal(j)=sens_normal(j)*(1-diagtest(d,t).sens(4)-diagtest(d,t).sens(5)); 
        else
            sens_normal(j)=sens_normal(j);
        end
    end
    for k=1:length(Above_norm)
        t=Above_norm(k);
        %sens_normal(j)=sens_normal(j)*(1-aboveSens(d,t)); 
        if ~isempty(diagtest(d,t).sens)
            sens_normal(j)=sens_normal(j)*(1-diagtest(d,t).sens(1)-diagtest(d,t).sens(2));
        else
            sens_normal(j)=sens_normal(j); 
        end
    end

    %continuous tests (above and below):
    for k=1:length(Below_b)  
        t=Below_b(k);
        %sens_medium(j,k,1)=1-belowSens(d,t);  %P(t= normal|d)
        %sens_medium(j,k,2)=1-ebelowSens(d,t);  %P(t>= below|d)
        if ~isempty(diagtest(d,t).sens)
            sens_medium(j,k,1)=1-diagtest(d,t).sens(4)-diagtest(d,t).sens(5);  %P(t= normal|d)
            sens_medium(j,k,2)=1-diagtest(d,t).sens(5);  %P(t>= below|d)
        else
            sens_medium(j,k,1)=1;  %P(t= normal|d)
            sens_medium(j,k,2)=1;  %P(t>= below|d)
        end
    end 
    for k=1:length(Above_a) 
        t=Above_a(k);
        %sens_medium(j,k+length(Below_b),1)=1-aboveSens(d,t);  %P(t= normal|d)
        %sens_medium(j,k+length(Below_b),2)=1-eaboveSens(d,t);  %P(t<= above|d)
        if ~isempty(diagtest(d,t).sens)
            sens_medium(j,k+length(Below_b),1)=1-diagtest(d,t).sens(1)-diagtest(d,t).sens(2);  %P(t= normal|d)
            sens_medium(j,k+length(Below_b),2)=1-diagtest(d,t).sens(1);  %P(t<= above|d)
        else
            sens_medium(j,k+length(Below_b),1)=1;  %P(t= normal|d)
            sens_medium(j,k+length(Below_b),2)=1;  %P(t<= above|d)
        end                
    end

     %continous tests (eabove and ebelow) + binary tests (true):
    for k=1:length(Binary_true)  
        t=Binary_true(k);
        sens_extreme(j,k)=1-restSens(d,t);  %P(t= normal|d)
    end
    for k=1:length(Below_eb) %noisyMAX: below > ebelow 
        t=Below_eb(k);
        %sens_extreme(j,k+length(Binary_true))=1-ebelowSens(d,t); %P(t>= below|d)
        if ~isempty(diagtest(d,t).sens)
            sens_extreme(j,k+length(Binary_true))=1-diagtest(d,t).sens(5); %P(t>= below|d)
        else
            sens_extreme(j,k+length(Binary_true))=1; %P(t>= below|d)
        end                
    end %for k
    for k=1:length(Above_ea) %noisyMAX: above < eabove 
        t=Above_ea(k);
        %sens_extreme(j,k+length(Binary_true)+length(Below_eb))=1-eaboveSens(d,t); %P(t<= above|d)
        if ~isempty(diagtest(d,t).sens)
            sens_extreme(j,k+length(Binary_true)+length(Below_eb))=1-diagtest(d,t).sens(1); %P(t<= above|d)
        else
            sens_extreme(j,k+length(Binary_true)+length(Below_eb))=1; %P(t<= above|d)
        end
    end %for k

end %for j 

%4. deleting from diagset diagnoses which do not allow normal valuue of
%negative tests:
index=find(sens_normal==0);
diagset(index)=[]; 
patientprev(index,:)=[];
sens_normal(index)=[];
sens_extreme(index,:)=[];
sens_medium(index,:,:)=[];

%5. implementing quickscore algorithm:
C1=de2bi([0:2^length(Medium)-1]); 
C2=de2bi([0:2^length(Extreme)-1]);
P_findings=0;  P_joint=zeros(length(diagset),3);   

dt_array = zeros(10,1); 
t1 = tic; 
if quickscore_run 
    for k1=1:size(C1,1)
    
        if size(C1,1)~=1
            sign1=(-1)^(sum(1-C1(k1,:)));
            S1=(prod(sens_medium(:,:,2).^repmat(C1(k1,:),length(diagset),1),2)).*(prod(sens_medium(:,:,1).^repmat((1-C1(k1,:)),length(diagset),1),2));
        else
            sign1=1;    
            S1=ones(length(diagset),1); 
        end
    
        for k2=1:size(C2,1)
            if size(C2,1)~=1
                if sum(C2(k2,:))==1
                    dt_array(log2(k2-1)+1) = toc(t1); 
                end
                sign2=(-1)^(sum(C2(k2,:)));
                S2=prod(sens_extreme.^repmat(C2(k2,:),length(diagset),1),2);
            else
                sign2=1;    
                S2=ones(length(diagset),1);
            end
    
            Term=sign1*sign2*prod(patientprev(:,1).*S1.*S2.*sens_normal+patientprev(:,2)); 
            P_findings=P_findings+Term; 
            
    %         display([C1(k1,:),C2(k2,:)])
    %         fprintf('%.16f ',S1(1:5)')
    %         fprintf('\n')
    %         fprintf('%.16f ',Term)
            for j=1:length(diagset)
               patientprev1=patientprev;
               patientprev1(j,2)=0; %prevalences to calculate P(dj=1,t1...tm)
               patientprev0=patientprev;
               patientprev0(j,1)=0; %prevalences to calculate P(dj=0,t1...tm) P(dj=1,t1...tm)
               Term1=sign1*sign2*prod(patientprev1(:,1).*S1.*S2.*sens_normal+patientprev1(:,2));
               Term2=sign1*sign2*prod(patientprev0(:,1).*S1.*S2.*sens_normal+patientprev0(:,2));
               P_joint(j,1)=P_joint(j,1)+Term1;
               P_joint(j,2)=P_joint(j,2)+Term2;
               if ismember(j,[1,2,3,4,5])
    %                fprintf('%.16f ',Term1)
               end
            end  %for j
    %         fprintf('\n\n\n')
        end %for k2
    end %for k1
end
dt_array(end) = toc(t1); 
dt_array = dt_array(2:end); 

P_joint(:,3)=P_joint(:,2)+P_joint(:,1)-P_findings;
if P_findings<10^-13
    fprintf('ERROR. Precision problem\n');
    fprintf(fid,'ERROR. Precision problem. P_findings= %f\n', P_findings);
end

pdiag=P_joint(:,1)/P_findings;  %P(di=1|t1...tm)=1 - P(di=0,t1...tm)/P(t1...tm)

[pdiag_sorted, index]=sort(pdiag,'descend');
diffdiag=diagset(index);
% P_joint=P_joint(index,:);

%6. Printing the results 
positivefindings=[Medium,Extreme];
negativefindings=setdiff(tests_1,positivefindings);
tests_2=patient(i).ext;
if ~isempty(tests_2)
    fprintf(fid,'%s\n','');
    fprintf(fid,'CONDITIONS:\n');
    for k=1:length(tests_2)
        fprintf(fid,'%-2d  t%-5d %-70s %s \n',k,tests_2(k),test(tests_2(k)).name,patient(i).evalue{k});
    end
end
if ~isempty(negativefindings)
    fprintf(fid,'%s\n','');
    fprintf(fid,'NEGATIVE TESTS:\n');
    for k=1:length(negativefindings)
        t=negativefindings(k);
        index=find(patient(i).test==t);
        value=patient(i).tvalue{index};
        if strcmp('normal',value)||strcmp('false',value)
            fprintf(fid,'%-2d  t%-5d %-70s %s\n',k,t,test(t).name,value);
        else
            fprintf(fid,'%-2d  %-70s %s\n',k,test(t).name,'ERROR:Inlavid value');
        end
    end
end
 if ~isempty(positivefindings)
    fprintf(fid,'%s\n','');
    fprintf(fid,'POSITIVE SYMPTOMS AND TESTS:\n');
    for k=1:length(positivefindings)
        t=positivefindings(k);
        index=find(patient(i).test==t);
        fprintf(fid,'%-2d  t%-5d %-70s %s\n',k,t,test(t).name,patient(i).tvalue{index});
    end
 end

Test_explain2=zeros(1,length(diffdiag));
Test_explain=cell(1,length(diffdiag));
for j=1:length(diffdiag)
    d=diffdiag(j);
    Test_explain{j}=zeros(1,length(positivefindings));
    for k=1:length(positivefindings)
        t=positivefindings(k);
        value=patient(i).tvalue{find(patient(i).test==t,1)};
        if ~isempty(diagtest(d,t).sens)
            if strmatch('below',value)
                if diagtest(d,t).sens(4)>0
                    Test_explain{j}(k)=diagtest(d,t).sens(4);
                end
            elseif strmatch('ebelow',value)
                if diagtest(d,t).sens(5)>0
                    Test_explain{j}(k)=diagtest(d,t).sens(5);
                end
            elseif strmatch('above',value)
                if diagtest(d,t).sens(2)>0
                    Test_explain{j}(k)=diagtest(d,t).sens(2);
                end
            elseif strmatch('eabove',value)
                if diagtest(d,t).sens(1)>0
                    Test_explain{j}(k)=diagtest(d,t).sens(1);
                end
            elseif restSens(d,t)>0
                Test_explain{j}(k)=restSens(d,t);
            end
        end
    end  
    Test_explain{j}=round(Test_explain{j},2);
    Test_explain2(j)=nnz(Test_explain{j});
end

Test_connect=zeros(length(positivefindings),1);
for j=1:length(positivefindings)
    t=positivefindings(j);
    value=patient(i).tvalue{find(patient(i).test==t,1)};
    for k=1:length(diffdiag)
        d=diffdiag(k);
        if ~isempty(diagtest(d,t).sens)
            if strmatch('below',value)
                if diagtest(d,t).sens(4)>0
                    Test_connect(j)=Test_connect(j)+1;
                end
            elseif strmatch('ebelow',value)
                if diagtest(d,t).sens(5)>0
                    Test_connect(j)=Test_connect(j)+1;
                end
            elseif strmatch('above',value)
                if diagtest(d,t).sens(2)>0
                    Test_connect(j)=Test_connect(j)+1;
                end
            elseif strmatch('eabove',value)
                if diagtest(d,t).sens(1)>0
                    Test_connect(j)=Test_connect(j)+1;
                end
            elseif restSens(d,t)>0
                Test_connect(j)=Test_connect(j)+1;
            end
        end
    end
end
formatSpec = '%6.0f'; 
Test_conn=num2str(Test_connect,formatSpec);
Test_conn=[repmat('[',length(positivefindings),1),Test_conn,repmat(']',length(positivefindings),1)];
Test_conn=[Test_conn,repmat(' ',length(positivefindings),6-size(Test_conn,2))];
Test_conn=reshape(Test_conn',1,[]);

Test_numb=num2str([1:length(positivefindings)],formatSpec);

fprintf(fid,'%s\n','');
fprintf(fid,'%-97s %s\n','KNOWN DIAGNOSIS','Explained positive tests:');
fprintf(fid,'Probab.   Preval.   Ndiag  Diagnosis %-62s %s\n','',Test_numb);
formatSpec = '%6.2f';
for j=1:length(patient(i).diag)
    index=find(diffdiag==patient(i).diag(j));
    Test_expl=num2str(Test_explain{index}',formatSpec);
    Test_expl(ismember(Test_expl,'0.00','rows'),:)=repmat('____',nnz(ismember(Test_expl,'0.00','rows')),1);
    Test_expl(:,end+1)=' '; Test_expl(:,end+1)=' ';
    Test_expl=reshape(Test_expl',1,[]);
    fprintf(fid,'%f  %-6f  d%-5d %-70s %s\n',pdiag_sorted(index),prev(diffdiag(index)),diffdiag(index),string(diagn(diffdiag(index))),Test_expl);
end

fprintf(fid,'%s\n','');
fprintf(fid,'%-97s %s\n','DIFFERENTIAL DIAGNOSIS (LARGEST PROBABILITIES)','Explained positive tests:');
fprintf(fid,'Probab.   Preval.   Ndiag  Diagnosis %-62s %s\n','',Test_numb);
fprintf(fid,'%-97s %s\n','Number of diagnoses connected to positive tests', Test_conn);
formatSpec = '%6.2f'; 
j=0; j1=0;
while j1<min(10,length(diagset))
    j=j+1;
    if ~ismember(diffdiag(j),patient(i).diag)
        Test_expl=num2str(Test_explain{j}',formatSpec);
        Test_expl(ismember(Test_expl,'0.00','rows'),:)=repmat('____',nnz(ismember(Test_expl,'0.00','rows')),1);
        Test_expl(:,end+1)=' '; Test_expl(:,end+1)=' ';
        Test_expl=reshape(Test_expl',1,[]);
        fprintf(fid,'%f  %-6f  d%-5d %-70s %s\n',pdiag_sorted(j),prev(diffdiag(j)),diffdiag(j),string(diagn(diffdiag(j))),Test_expl); 
    j1=j1+1;
    end
end


[~, index]=sort(Test_explain2,'descend');
pdiag_sorted2=pdiag_sorted(index);
diffdiag2=diffdiag(index);
Test_explain=Test_explain(index);

fprintf(fid,'%s\n','');
fprintf(fid,'%-97s %s\n','DIFFERENTIAL DIAGNOSIS (BEST POSITIVE TESTS EXPLANATION)','Explained positive tests:');
fprintf(fid,'Probab.   Preval.   Ndiag  Diagnosis %-62s %s\n','',Test_numb);

for j=1:min(10,length(diagset)) 
    Test_expl=num2str(Test_explain{j}',formatSpec);
    Test_expl(ismember(Test_expl,'0.00','rows'),:)=repmat('____',nnz(ismember(Test_expl,'0.00','rows')),1);
    Test_expl(:,end+1)=' '; Test_expl(:,end+1)=' ';
    Test_expl=reshape(Test_expl',1,[]);
    fprintf(fid,'%f  %-6f  d%-5d %-70s %s\n',pdiag_sorted2(j),prev(diffdiag2(j)),diffdiag2(j),string(diagn(diffdiag2(j))),Test_expl); 
end

fprintf(fid,'%s\n','');
fprintf(fid,'%-97s %s\n','CORRECT DIAGNOSIS (GIVEN)','Explained positive tests:');
fprintf(fid,'Probab.   Preval.   Ndiag  Diagnosis %-62s %s\n','',Test_numb);
correct=patient(i).cdiag;
correct(find(correct==0))=[];
for j=1:length(correct)
    d=correct(j);
    Test_explain=zeros(1,length(positivefindings));
    for k=1:length(positivefindings)
        t=positivefindings(k);
        value=patient(i).tvalue{find(patient(i).test==t,1)};
        if ~isempty(diagtest(d,t).sens)
            if strmatch('below',value)
                if diagtest(d,t).sens(4)>0
                    Test_explain(k)=diagtest(d,t).sens(4);
                end
            elseif strmatch('ebelow',value)
                if diagtest(d,t).sens(5)>0
                    Test_explain(k)=diagtest(d,t).sens(5);
                end
            elseif strmatch('above',value)
                if diagtest(d,t).sens(2)>0
                    Test_explain(k)=diagtest(d,t).sens(2);
                end
            elseif strmatch('eabove',value)
                if diagtest(d,t).sens(1)>0
                    Test_explain(k)=diagtest(d,t).sens(1);
                end
            elseif restSens(d,t)>0
                Test_explain(k)=restSens(d,t);
            end
        end
    end  
    Test_explain=round(Test_explain,2);

    Test_expl=num2str(Test_explain',formatSpec);
    Test_expl(ismember(Test_expl,'0.00','rows'),:)=repmat('____',nnz(ismember(Test_expl,'0.00','rows')),1);
    Test_expl(:,end+1)=' '; Test_expl(:,end+1)=' ';
    Test_expl=reshape(Test_expl',1,[]);

    index=find(diffdiag==correct(j));
    if index>0
        fprintf(fid,'%f  %-6f  d%-5d %-70s %s\n',pdiag_sorted(index),prev(correct(j)),correct(j),string(diagn(correct(j))),Test_expl); 
    else
        fprintf(fid,'--------  %-6f  d%-5d %-70s %s\n',prev(correct(j)),correct(j),string(diagn(correct(j))),Test_expl);
    end
end
    fprintf(fid,'%s\n',''); 

end %if ~isempty(patient(i).age)

timeval=toc;
fprintf(fid,'%s %f %s\n','Time= ',timeval,' seconds.');
fprintf(fid,'%s\n','');fprintf(fid,'%s\n',''); fprintf(fid,'%s\n',''); 




fclose(fid);

