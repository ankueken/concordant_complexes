% data for supplementary table S1

clear

files_irrev=dir('Results_concordant_original_irreversibility_considered\*.mat');
files_obj=dir('Results_concordant_original_objective\*.mat');

Data_irrev=[];Data_obj=[];Data_in_concordant_irrev=[];Data_in_concordant_irrev_rel=[];
Data_irrev_rel=[];Data_obj_rel=[];Data_in_concordant_obj=[];Data_in_concordant_obj_rel=[];

for f=1:length(files_irrev)
    
    R_irrev=load(strcat(files_irrev(f).folder,'/',files_irrev(f).name));
    R_obj=load(strcat(files_obj(f).folder,'/',files_obj(f).name));
       
    CC_triu = triu(R_irrev.CC);
    total_pairs=sum(sum((triu(ones(size(R_irrev.CC)),1))));

    balanced_concordant = sum(sum(CC_triu==1));
    trivial_concordant = sum(sum(CC_triu==-1));
    higher_order_concordant_pairs = sum(sum(CC_triu==2));
    total=balanced_concordant+trivial_concordant+higher_order_concordant_pairs;
    positive_ratio = length(intersect(find(R_irrev.coupling_pairs(:,3)>0),find(R_irrev.coupling_pairs(:,3)<1000)));
    [in_concordant_irrev,~] = find(CC_triu~=0);
    
    Data_irrev=[Data_irrev;balanced_concordant trivial_concordant higher_order_concordant_pairs positive_ratio total];
    Data_irrev_rel=[Data_irrev_rel; balanced_concordant/total_pairs trivial_concordant/total_pairs higher_order_concordant_pairs/total_pairs positive_ratio/total_pairs total/total_pairs];
    Data_in_concordant_irrev=[Data_in_concordant_irrev; length(unique(in_concordant_irrev))];
    Data_in_concordant_irrev_rel=[Data_in_concordant_irrev_rel; length(unique(in_concordant_irrev))/size(CC_triu,1)];
    
    CC_triu = triu(R_obj.CC);
    total_pairs=sum(sum((triu(ones(size(R_obj.CC)),1))));

    balanced_concordant = sum(sum(CC_triu==1));
    trivial_concordant = sum(sum(CC_triu==-1));
    higher_order_concordant_pairs = sum(sum(CC_triu==2));
    total=balanced_concordant+trivial_concordant+higher_order_concordant_pairs;
     positive_ratio = length(intersect(find(R_obj.coupling_pairs(:,3)>0),find(R_obj.coupling_pairs(:,3)<1000)));
    [in_concordant_obj,~] = find(CC_triu~=0);
        
    Data_obj=[Data_obj;balanced_concordant trivial_concordant higher_order_concordant_pairs positive_ratio total];
    Data_obj_rel=[Data_obj_rel; balanced_concordant/total_pairs trivial_concordant/total_pairs higher_order_concordant_pairs/total_pairs positive_ratio/total_pairs total/total_pairs];
    Data_in_concordant_obj=[Data_in_concordant_obj; length(unique(in_concordant_obj))];
    Data_in_concordant_obj_rel=[Data_in_concordant_obj_rel; length(unique(in_concordant_obj))/size(CC_triu,1)];
    Number_complexes_and_pairs(f,1:2) = [size(CC_triu,1) total_pairs];
end

Data_irrev_rel=Data_irrev_rel*100;
Data_obj_rel=Data_obj_rel*100;
Data_in_concordant_irrev_rel=Data_in_concordant_irrev_rel*100;
Data_in_concordant_obj_rel=Data_in_concordant_obj_rel*100;

% model_names={'\itA. niger iMA871';
%     '\itA. thaliana AraCore';
%     '\itC. reinhardtii iCre1355';
%     '\itE. coli iJO1366';
%     '\itM. acetivorans iMB745';
%     '\itM. barkeri iAF692';
%     '\itM. musculus';
%     '\itM. tuberculosis iNJ661m';
%     '\itN. pharaonis';
%     '\itP. putida iJN746';
%     '\itT. maritima iLJ478';
%     '\itS. cerevisiae Yeast8'};
% 
% % Figure 2
% subplot(2,1,1)
% bar(1:3:36,Data_irrev_rel(:,[2 3 1]),0.3,'stacked')
% hold on
% bar(2:3:36,Data_obj_rel(:,[2 3 1]),0.3,'stacked')
% set(gca,'YScale','Log')
% ylabel('Concordant complex pairs (%)')
% set(gca,'XTick',1.5:3:36,'XTickLabels',[],'XTickLabelRotation',45)
% ylim([0 100])
% 
% subplot(2,1,2)
% bar(1:3:36,Data_in_concordant_irrev_rel,0.3,'stacked')
% hold on
% bar(2:3:36,Data_in_concordant_obj_rel,0.3,'stacked')
% ylabel('Complexes in concordance relation (%)')
% set(gca,'XTick',1.5:3:36,'XTickLabels',model_names,'XTickLabelRotation',45)
% ylim([0 100])

