% calculate effective degree

files_irrev=dir('Results_concordant_original_irreversibility_considered\*.mat');
files_obj=dir('Results_concordant_original_objective\*.mat');

mkdir('effective_degreeseq_irreversible')
mkdir('effective_degreeseq_objective')

effdeg_irrev=[];effdeg_obj=[];

for f=1:length(files_irrev)

    R_irrev=load(strcat(files_irrev(f).folder,'/',files_irrev(f).name));

    complex_class_matrix=zeros(size(R_irrev.CC,1),length(R_irrev.class_with_balanced));

    for i=1:size(R_irrev.CC,1)
        for j=1:length(R_irrev.class_with_balanced)
            complex_class_matrix(i,j)=any(ismember(R_irrev.class_with_balanced{j},i));        
        end
    end

    for s=1:length(R_irrev.species_degree_complexes)
        effective_degree(s,1)=sum(sum(complex_class_matrix(find(R_irrev.Results_balanced.MODEL_r{1}.Y(s,:)~=0),:))~=0);
    end
    
    [H,b]=hist(effective_degree,1:max(effective_degree));
    b(H==0)=[];
    H(H==0)=[];
    T=table(b',H','VariableNames',{'xvalue','counts'});
    effdeg_irrev=[effdeg_irrev; mean(effective_degree) min(effective_degree) max(effective_degree)];
    writetable(T,strcat('effective_degreeseq_irreversible/',files_irrev(f).name(1:end-4),'.txt'),'Delimiter',',')
    clear T R_irrev effective_degree 
    
    %%
    R_obj=load(strcat(files_obj(f).folder,'/',files_obj(f).name));

    complex_class_matrix=zeros(size(R_obj.CC,1),length(R_obj.class_with_balanced));

    for i=1:size(R_obj.CC,1)
        for j=1:length(R_obj.class_with_balanced)
            complex_class_matrix(i,j)=any(ismember(R_obj.class_with_balanced{j},i));        
        end
    end

    for s=1:length(R_obj.species_degree_complexes)
        effective_degree(s,1)=sum(sum(complex_class_matrix(find(R_obj.Results_balanced.MODEL_o{1}.Y(s,:)~=0),:))~=0);
    end
    
    [H,b]=hist(effective_degree,1:max(effective_degree));
    b(H==0)=[];
    H(H==0)=[];
    T=table(b',H','VariableNames',{'xvalue','counts'});
    effdeg_obj=[effdeg_obj; mean(effective_degree) min(effective_degree) max(effective_degree)];
    writetable(T,strcat('effective_degreeseq_objective/',files_obj(f).name(1:end-4),'.txt'),'Delimiter',',')
    clear T R_obj effective_degree 
end