% analyse statistics for random variants in comparison to E. coli model

clear
files=dir('Results_concordant_randomization_irreversibility_considered\*.mat');

for f=1:length(files)
    
    R=load(strcat(files(f).folder,'/',files(f).name));
    
% reducability index
    RI_no_balanced(f,1) = 1-(length(R.class)/(size(R.CC,1)-length(R.B)));
    
% number of concordance modules
    number_modules(f,1) = length(R.class);
   
% the size of the largest concordance module
    size_largest_module(f,1) = max(cellfun(@length,R.class));

% the mean concordance module size
    mean_size_module(f,1) = mean(cellfun(@length,R.class));

% --------------------------------------------------   
    
     for s=1:size(R.Results_balanced.MODEL_r{1}.Y,1)
        complexes_of_species = find(R.Results_balanced.MODEL_r{1}.Y(s,:)~=0);
        for c=1:length(R.class)
            species_module_matrix{f}(s,c)=~isempty(intersect(complexes_of_species,R.class{c})); 
        end
     end
     
% the mean number of concordance modules in which a metabolite participates
    mean_number_modules_per_metabolite(f,1) = mean(sum(species_module_matrix{f},2));

% the maximum number of concordance modules in which metabolites participate
    max_number_modules_per_metabolite(f,1) = max(sum(species_module_matrix{f},2));

end

RO = load('Results_concordant_original_irreversibility_considered\Ecoli2011-iJO1366.mat');

% reducability index
    RIO_no_balanced = 1-(length(RO.class)/(size(RO.CC,1)-length(RO.B)));

% number of concordance modules
    RO_number_modules = length(RO.class);
   
% the size of the largest concordance module
    RO_size_largest_module = max(cellfun(@length,RO.class));

% the mean concordance module size
    RO_mean_size_module = mean(cellfun(@length,RO.class));

% --------------------------------------------------   
    
     for s=1:size(RO.Results_balanced.MODEL_r{1}.Y,1)
        complexes_of_species = find(RO.Results_balanced.MODEL_r{1}.Y(s,:)~=0);
        for c=1:length(RO.class)
            RO_species_module_matrix(s,c)=~isempty(intersect(complexes_of_species,RO.class{c})); 
        end
     end
    
% the mean number of concordance modules in which a metabolite participates
    RO_mean_number_modules_per_metabolite = mean(sum(RO_species_module_matrix,2));

% the maximum number of concordance modules in which metabolites participate
    RO_max_number_modules_per_metabolite = max(sum(RO_species_module_matrix,2));

% z-score 
    
% number of concordance modules
    z_number_modules = (RO_number_modules-mean(number_modules))/std(number_modules);
    p_number_modules = 2*normcdf(-abs(z_number_modules));

% reducability index
    z_RI = (RIO_no_balanced-mean(RI_no_balanced))/std(RI_no_balanced);
    p_RI = 2*normcdf(-abs(z_RI)); 
     
% the size of the largest concordance module
    z_size_largest_module = (RO_size_largest_module-mean(size_largest_module))/std(size_largest_module);
    p_size_largest_module = 2*normcdf(-abs(z_size_largest_module)); 
    
% the mean concordance module size
    z_mean_size_module = (RO_mean_size_module-mean(mean_size_module))/std(mean_size_module);
    p_mean_size_module = 2*normcdf(-abs(z_mean_size_module)); 
    
% the mean number of concordance modules in which a metabolite participates
    z_mean_number_modules_per_metabolite = (RO_mean_number_modules_per_metabolite-mean(mean_number_modules_per_metabolite))/std(mean_number_modules_per_metabolite);
    p_mean_number_modules_per_metabolite = 2*normcdf(-abs(z_mean_number_modules_per_metabolite)); 
    
% the maximum number of concordance modules in which metabolites participate
    z_max_number_modules_per_metabolite = (RO_max_number_modules_per_metabolite-mean(max_number_modules_per_metabolite))/std(max_number_modules_per_metabolite);
    p_max_number_modules_per_metabolite = 2*normcdf(-abs(z_max_number_modules_per_metabolite)); 
    
figure
    subplot(2,3,1)
% number of modules
    histogram(number_modules,10,'Normalization','probability')
    hold on
    plot([RO_number_modules RO_number_modules],[0 1])
    ylabel('Probability')
    xlabel('Number of concordance modules')
    legend('randomized','iJO1366')
    
    subplot(2,3,2)
% the size of the largest concordance module
    histogram(size_largest_module,10,'Normalization','probability')
    hold on
    plot([RO_size_largest_module RO_size_largest_module],[0 1])
    ylabel('Probability')
    xlabel('Size of largest concordance module')
    
    subplot(2,3,3)
% the mean concordance module size
    histogram(mean_size_module,10,'Normalization','probability')
    hold on
    plot([RO_mean_size_module RO_mean_size_module],[0 1])
    ylabel('Probability')
    xlabel('Average size of concordance module')
    
    subplot(2,3,4)
% the mean number of concordance modules in which a metabolite participates
    histogram(mean_number_modules_per_metabolite,10,'Normalization','probability')
    hold on
    plot([RO_mean_number_modules_per_metabolite RO_mean_number_modules_per_metabolite],[0 1])
    ylabel('Probability')
    xlabel('Average number of modules per metabolite')

    subplot(2,3,5)
% the maximum number of concordance modules in which metabolites participate
    histogram(max_number_modules_per_metabolite,10,'Normalization','probability')
    hold on
    plot([RO_max_number_modules_per_metabolite RO_max_number_modules_per_metabolite],[0 1])
    ylabel('Probability')
    xlabel('Maximum number of modules per metabolite')

    subplot(2,3,6)
% reducability index
    histogram(RI_no_balanced,10,'Normalization','probability')
    hold on
    plot([RIO_no_balanced RIO_no_balanced],[0 1])
    ylabel('Probability')
    xlabel('Reducability index')
    
    