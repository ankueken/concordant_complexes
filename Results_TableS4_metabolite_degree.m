% supplementary table S4
clear

model_names={'A. niger iMA871';
    'A. thaliana AraCore';
    'C. reinhardtii iCre1355';
    'E. coli iJO1366';
    'M. acetivorans iMB745';
    'M. barkeri iAF692';
    'M. musculus';
    'M. tuberculosis iNJ661m';
    'N. pharaonis';
    'P. putida iJN746';
    'T. maritima iLJ478';
    'S. cerevisiae Yeast8'};

files=dir('Results_concordant_original_irreversibility_considered\*test.mat');
T=table();
for f=1:length(files)
    R=load(strcat(files(f).folder,'/',files(f).name));

    Tn=table(R.Results_balanced.MODEL_r{1}.mets,  R.Results_balanced.MODEL_r{1}.metNames, ...
        repmat(model_names(f),size(R.species_degree_complexes)), R.species_degree_complexes,...
        'VariableNames',{'metabolite abbreviation','metabolite full name','organism','degree'});
    T=[T;Tn];
end

writetable(T,'SupplementaryTableS4.xlsx')
    
T_ORDERED = sortrows(T,4,'descend');

