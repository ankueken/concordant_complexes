clear

files_irrev=dir('Results_concordant_original_irreversibility_considered\*.mat');
files_obj=dir('Results_concordant_original_objective\*.mat');

for f=1:length(files_irrev)
    
    R_irrev=load(strcat(files_irrev(f).folder,'/',files_irrev(f).name));
    R_obj=load(strcat(files_obj(f).folder,'/',files_obj(f).name));
    
    RI_irrev(f,1) = 1-(length(R_irrev.class_with_balanced)/size(R_irrev.CC,1));
    RI_irrev_no_balanced(f,1) = 1-(length(R_irrev.class)/(size(R_irrev.CC,1)-length(R_irrev.B)));
    
    RI_obj(f,1) = 1-(length(R_obj.class_with_balanced)/size(R_obj.CC,1));
    RI_obj_no_balanced(f,1) = 1-(length(R_obj.class)/(size(R_obj.CC,1)-length(R_obj.B)));
end

model_names={'\itA. niger iMA871';
    '\itA. thaliana AraCore';
    '\itC. reinhardtii iCre1355';
    '\itE. coli iJO1366';
    '\itM. acetivorans iMB745';
    '\itM. barkeri iAF692';
    '\itM. musculus';
    '\itM. tuberculosis iNJ661m';
    '\itN. pharaonis';
    '\itP. putida iJN746';
    '\itT. maritima iLJ478';
    '\itS. cerevisiae Yeast8'};

% figure S5
[~,idx] = sort(RI_irrev);
bar(1:3:36,RI_irrev(idx),0.3)
hold on
bar(2:3:36,RI_obj(idx),0.3)
ylabel('Reducibility index')
set(gca,'XTick',1.5:3:36,'XTickLabels',model_names(idx),'XTickLabelRotation',45)
legend('(i) irreversibility constraints','(ii) optimality constraint','Orientation','horizontal')
legend boxoff

% figure 5
[~,idx] = sort(RI_irrev_no_balanced);
bar(1:3:36,RI_irrev_no_balanced(idx),0.3)
hold on
bar(2:3:36,RI_obj_no_balanced(idx),0.3)
ylabel('Reducibility index')
set(gca,'XTick',1.5:3:36,'XTickLabels',model_names(idx),'XTickLabelRotation',45)
legend('(i) irreversibility constraints','(ii) optimality constraint','Orientation','horizontal')
legend boxoff
