% subject demographics t test

load('/Users/carlynpattersongentile/Documents/MATLAB/Carlyn/vepMELA_subjectInfo.mat')

MwA=[3 6 8 9 10 14 17 18]';

HAf=[1 2 4 5 7 11 12 13 15 16]';

temp=table2array(scoreTable(MwA,{'Age'}));
temp2=cellstr(temp);
temp3=cell2mat(temp2);
temp4=temp3(:,1:2);
MwA_age=str2num(temp4);
clear temp*

temp=table2array(scoreTable(HAf,{'Age'}));
temp2=cellstr(temp);
temp3=cell2mat(temp2);
temp4=temp3(:,1:2);
HAf_age=str2num(temp4);
clear temp*

[h_age,p_age]=ttest2(MwA_age,HAf_age);

temp=table2array(scoreTable(MwA,{'Sex'}));
MwA_sex=zeros(1,length(temp));
for x=1:length(temp)
    if temp(x)=='Female'
        MwA_sex(:,x)=1;
    end
end
clear temp* x

temp=table2array(scoreTable(HAf,{'Sex'}));
HAf_sex=zeros(1,length(temp));
for x=1:length(temp)
    if temp(x)=='Female'
        HAf_sex(:,x)=1;
    end
end
clear temp* x

[h_sex,p_sex]=ttest2(MwA_sex,HAf_sex);


[h_VDS,p_VDS]=ttest2(table2array(scoreTable(MwA,{'VDS'})),table2array(scoreTable(HAf,{'VDS'})));