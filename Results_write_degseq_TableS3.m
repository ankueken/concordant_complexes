%% create files of metabolite degree, module size (Table S3)

files_irrev=dir('Results_concordant_original_irreversibility_considered\*.mat');
files_obj=dir('Results_concordant_original_objective\*.mat');

mkdir('degseq')
mkdir('moduleseq_irreversible')
mkdir('moduleseq_objective')
mkdir('moduleseq_no_balanced_irreversible')
mkdir('moduleseq_no_balanced_objective')

modseq_irrev=[]; modseq_wb_irrev=[];degseq=[];
modseq_obj=[]; modseq_wb_obj=[];

for f=1:length(files_irrev)

   R_irrev=load(strcat(files_irrev(f).folder,'/',files_irrev(f).name));
   [H,b]=hist(R_irrev.species_degree_complexes,1:max(R_irrev.species_degree_complexes));
   b(H==0)=[];
   H(H==0)=[];
   T=table(b',H','VariableNames',{'xvalue','counts'});
   degseq = [degseq;mean(R_irrev.species_degree_complexes) min(R_irrev.species_degree_complexes) max(R_irrev.species_degree_complexes)];
   writetable(T,strcat('degseq/',files_irrev(f).name(1:end-4),'.txt'),'Delimiter',',')
   
   [H,b]=hist(cellfun(@length,R_irrev.class_with_balanced),1:max(cellfun(@length,R_irrev.class_with_balanced)));
   b(H==0)=[];
   H(H==0)=[];
   T=table(b',H','VariableNames',{'xvalue','counts'});
   modseq_wb_irrev = [modseq_wb_irrev;length(R_irrev.class_with_balanced) mean(cellfun(@length,R_irrev.class_with_balanced)) min(cellfun(@length,R_irrev.class_with_balanced)) max(cellfun(@length,R_irrev.class_with_balanced))];
   writetable(T,strcat('moduleseq_irreversible/',files_irrev(f).name(1:end-4),'.txt'),'Delimiter',',')
   
   [H,b]=hist(cellfun(@length,R_irrev.class),1:max(cellfun(@length,R_irrev.class)));
   b(H==0)=[];
   H(H==0)=[];
   T=table(b',H','VariableNames',{'xvalue','counts'});
   modseq_irrev = [modseq_irrev;length(R_irrev.class) mean(cellfun(@length,R_irrev.class)) min(cellfun(@length,R_irrev.class)) max(cellfun(@length,R_irrev.class))];
   writetable(T,strcat('moduleseq_no_balanced_irreversible/',files_irrev(f).name(1:end-4),'.txt'),'Delimiter',',')
   
   clear T R_irrev
    
 %%   
   R_obj=load(strcat(files_obj(f).folder,'/',files_obj(f).name));
   
   
   [H,b]=hist(cellfun(@length,R_obj.class_with_balanced),1:max(cellfun(@length,R_obj.class_with_balanced)));
   b(H==0)=[];
   H(H==0)=[];
   T=table(b',H','VariableNames',{'xvalue','counts'});
   modseq_wb_obj = [modseq_wb_obj;length(R_obj.class_with_balanced) mean(cellfun(@length,R_obj.class_with_balanced)) min(cellfun(@length,R_obj.class_with_balanced)) max(cellfun(@length,R_obj.class_with_balanced))];
   writetable(T,strcat('moduleseq_objective/',files_obj(f).name(1:end-4),'.txt'),'Delimiter',',')
   
   [H,b]=hist(cellfun(@length,R_obj.class),1:max(cellfun(@length,R_obj.class)));
   b(H==0)=[];
   H(H==0)=[];
   T=table(b',H','VariableNames',{'xvalue','counts'});
   modseq_obj = [modseq_obj;length(R_obj.class) mean(cellfun(@length,R_obj.class)) min(cellfun(@length,R_obj.class)) max(cellfun(@length,R_obj.class))];
   writetable(T,strcat('moduleseq_no_balanced_objective/',files_obj(f).name(1:end-4),'.txt'),'Delimiter',',')
   
   clear T R_obj
   
end
