% script to classify concordant modules

clear
addpath('functions')
files=dir('Results_concordant_original_irreversibility_considered\*.mat');

hit_species_single_module=cell(12,1);

for f=1%:length(files)
    disp(f)
    
    clearvars -except f files module_class percentage_module_classes number_module_classes ...
    hit_species_single_module independent_modules number_pseudo_independent_modules number_independent_modules ...
    R1 R1_closed R2 R3 R4 R5 module_size percentage_pseudo_independent_modules size_closed percentage_independent_modules

    load(strcat(files(f).folder,'/',files(f).name)) 
    size_closed{f} = nan(length(class),1);
    
%%   classification of concordant modules 

%   take networks without balanced networks and check input/output for each
%   module in that network:

    Results_balanced.MODEL_r{2} = remove_balanced_complexes_any(Results_balanced.MODEL_r{1},B);
    
    for c=1:length(class)
         sub_A = Results_balanced.MODEL_r{2}.A(class{c},:); % subnetwork/module
         incoming = ~isempty(setdiff(find(all(sub_A>=0)), find(all(sub_A==0))));
         outgoing = ~isempty(setdiff(find(all(sub_A<=0)), find(all(sub_A==0))));
        
%   (1) source concordant modules, 
%   that have no input from any complex outside of the concordant module, 
%   but have output to other concordant modules

        if incoming == 0 && outgoing == 1
            module_class{f}(c,1) = 1;

%   (2) sink concordant modules, 
%   that have no output to any complex outside of the concordant module, 
%   but have some inputs from other concordant modules

        elseif incoming == 1 && outgoing == 0
            module_class{f}(c,1) = 2;

%   (3) intermediate concordant modules, 
%   that have input and output from complexes outside of the concordant module

        elseif incoming == 1 && outgoing == 1
            module_class{f}(c,1) = 3;

%   (4) closed concordant modules, 
%   that have no input or output from any complex outside of the concordant module

        elseif incoming == 0 && outgoing == 0
            module_class{f}(c,1) = 4;
            size_closed{f}(c,1) = length(class{c});
            
        else
            module_class{f}(c,1) = 0;
        end
    end   
    
    percentage_module_classes(1,f) = (length(find(module_class{f}==1))/length(class))*100; 
    percentage_module_classes(2,f) = (length(find(module_class{f}==2))/length(class))*100; 
    percentage_module_classes(3,f) = (length(find(module_class{f}==3))/length(class))*100; 
    percentage_module_classes(4,f) = (length(find(module_class{f}==4))/length(class))*100; 
    
    number_module_classes(1,f) = length(find(module_class{f}==1)); 
    number_module_classes(2,f) = length(find(module_class{f}==2)); 
    number_module_classes(3,f) = length(find(module_class{f}==3)); 
    number_module_classes(4,f) = length(find(module_class{f}==4)); 
    
%   An input-closed module, i.e. a source or closed module, that includes 
%   species that do not appear in any other module in the network will be 
%   called independent

    for s=1:size(Results_balanced.MODEL_r{2}.Y,1)
        complexes_of_species = find(Results_balanced.MODEL_r{2}.Y(s,:)~=0);
        for c=1:length(class)
            species_module_matrix{f}(s,c)=~isempty(intersect(complexes_of_species,class{c})); 
        end
    end
    
    number_module_species = sum(species_module_matrix{f},2);
    
    independent_modules{f} = zeros(length(class),1);
    
    % all species buffered ?
    for i=1:size(species_module_matrix,2)
        if all(species_degree_complexes(species_module_matrix{f}(:,i)~=0)>=quantile(species_degree_complexes,.9))
            independent_modules{f}(i,1) = 2;
        end
    end
    
    % all single module species ?
    for i=1:size(species_module_matrix,2)
        if all(number_module_species(species_module_matrix{f}(:,i)~=0)==1)
            independent_modules{f}(i,1) = 1;
        end
    end
    
    % single module species and buffered species ?
    for s=1:size(Results_balanced.MODEL_r{2}.Y,1)
        if sum(species_module_matrix{f}(s,:))==1
            modules_single_module_species = find(species_module_matrix{f}(s,:)~=0);
            species_single_species_module = find(species_module_matrix{f}(:,modules_single_module_species));
            
            if isempty(setdiff(species_single_species_module,s))
                independent_modules{f}(modules_single_module_species) = 1;
                
            elseif species_degree_complexes(setdiff(species_single_species_module,s))>=quantile(species_degree_complexes,.9)
                independent_modules{f}(modules_single_module_species) = 2;
                hit_species_single_module{f}(end+1,:)=[Results_balanced.MODEL_r{2}.mets(s),s,modules_single_module_species];
            end
        end
    end
    
    %% Results Table S5 and S6
    % size modules
    module_size{f}=cellfun(@length,class);
    % size of pseudo-independent modules
    R1{f}=module_size{f}(intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==2)))
    % size of pseudo-independent modules
    R1_closed(f,:)=[min(module_size{f}(find(module_class{f}==4)))...
        max(module_size{f}(find(module_class{f}==4)))...
        mean(module_size{f}(find(module_class{f}==4)))...
        median(module_size{f}(find(module_class{f}==4)))];
    % min/max/mean number of metabolites in a pseudo-independent module
    R2(f,:)=[min(sum(species_module_matrix{f}(:,intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==2)))))...
        max(sum(species_module_matrix{f}(:,intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==2)))))...
        mean(sum(species_module_matrix{f}(:,intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==2)))))];
    
    m=intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==2));
    
    for i=1:length(m)
        [~,rxns_in_module] = find(Results_balanced.MODEL_r{1}.A(class{m(i)},:)~=0);
        Results_balanced.MODEL_r{1}.subSystems(rxns_in_module)
        [rows,~]=find(species_module_matrix{f}(:,m(i)));
        
        % Table S6 - metabolites in pseudo-independent modules
        % results shown in table S6 is the combination of these itterative
        % results, here we loop over modules and do not store them
        R4=[Results_balanced.MODEL_r{1}.mets(rows) Results_balanced.MODEL_r{1}.metNames(rows)];
        R5=number_module_species(rows);
    end
    
    number_independent_modules(f,1) = length(intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==1)));
    number_pseudo_independent_modules(f,1) = length(intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==2)));

    percentage_independent_modules(f,1) = (length(intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==1)))/length(class))*100;
    percentage_pseudo_independent_modules(f,1) = (length(intersect(union(find(module_class{f}==4),find(module_class{f}==1)),find(independent_modules{f}==2)))/length(class))*100;

end

clearvars -except f files module_class percentage_module_classes number_module_classes ...
    hit_species_single_module independent_modules number_pseudo_independent_modules number_independent_modules ...
    R1 R1_closed R2 R3 R4 R5 module_size percentage_pseudo_independent_modules size_closed percentage_independent_modules

