%% Network prepreparation, find balanced complexes and concordant complex pairs
% irreversibility constraints considered

% To run the code please adapt Line 23 -> path to R (third argument)
% Results_balanced.MODEL_o{1}=clean_model(name,model,'objective','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"',1e-6);

clear

addpath('functions')

species_list={'A_niger_iMA871';'ArabidopsisCoreModel';'M_acetivorans_iMB745';
    'Ecoli2011-iJO1366';'M_musculus';'M_barkeri_iAF692';'T_maritima_iLJ478';
    'C_reinhardtii_iCre1355_auto';'M_tuberculosis_iNJ661m';'P_putida_iJN746';'N_pharaonis';'YeastGEM'};

for f_n=1:length(species_list)
    clc
    clearvars -except f_n species_list
    disp(f_n)
    name = species_list{f_n};
    load(strcat('GEM_mat/',name,'.mat'))
% 
%     %% create models for each scenario
    Results_balanced.MODEL_r{1}=clean_model(name,model,'irreversibility_considered','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"',1e-6);

    %% species degree

    species_degree_complexes = sum(Results_balanced.MODEL_r{1}.Y~=0,2);
    species_degree_reactions = sum(Results_balanced.MODEL_r{1}.S~=0,2);

    %% balanced complexes
    [B,group] = find_balanced_complexes(Results_balanced.MODEL_r{1},1e-6);

%% 
% FIND PAIRS OF CONCORDANT COMPLEXES
% 
% If for two complexes $c_i$ and $c_j$ it hold that $\frac{A_i v}{A_j v}=\gamma$ 
% at any steady state, then we call those two complexes to be concordant.
% 
% If complexes $c_i$ and $c_j$ share a degree 2 species, then we call them to 
% be trivially concordant.
% 
% _Each pair is counted only once we either store ratio_ $\frac{i}{j}$ _or_ 
% $\frac{j}{i}$, not both.
% 
%     %% 1 - mutually concordant balanced complexes
% 
    CC = zeros(length(Results_balanced.MODEL_r{1}.complexes));
    CC(B,B) = 1;
% 
%     %% -1 - trivial concordant complexes
%     % complexes including degree 2 species
% 
    idx = find(species_degree_complexes==2);
    for i=1:length(idx)
        tcc = find(Results_balanced.MODEL_r{1}.Y(idx(i),:)~=0);
        if CC(tcc,tcc)~=1
            CC(tcc,tcc) = -1;
            CC(logical(eye(size(CC)))) = 0;
        end
    end

    %% 2 - other couplings

    coupling_pairs = find_coupling_complexes(Results_balanced.MODEL_r{1},B,f_n,0,group);

    for i = 1:size(coupling_pairs,1)
        if CC(coupling_pairs(i,1),coupling_pairs(i,2)) == 0
            CC(coupling_pairs(i,1),coupling_pairs(i,2)) = 2;
            CC(coupling_pairs(i,2),coupling_pairs(i,1)) = 2;
        end
    end

%% 
% *Reaction pattern around concordant complexes*

    reactions_concordant_complexes = array2table(nan(length(coupling_pairs),4));
    reactions_concordant_complexes.Properties.VariableNames = {'Ci_in','Ci_out','Cj_in','Cj_out'};
    
    for i=1:size(coupling_pairs,1)
        reactions_concordant_complexes.Ci_in(i,1) = sum(Results_balanced.MODEL_r{1}.A(coupling_pairs(i,1),:)>0);
        reactions_concordant_complexes.Ci_out(i,1) = sum(Results_balanced.MODEL_r{1}.A(coupling_pairs(i,1),:)<0);
        reactions_concordant_complexes.Cj_in(i,1) = sum(Results_balanced.MODEL_r{1}.A(coupling_pairs(i,2),:)>0);
        reactions_concordant_complexes.Cj_out(i,1) = sum(Results_balanced.MODEL_r{1}.A(coupling_pairs(i,2),:)<0);
    end
    
%% 
% *Concordance modules*
    
    % *Group mutually concordant complexes without balanced*

   [CP1(:,1),CP1(:,2)] = find(CC==-1);
   [CP2(:,1),CP2(:,2)] = find(CC==2);
   CP = [CP1; CP2];
   unclassified = unique(CP);

    class=[];
    while ~isempty(unclassified)
        i = unclassified(1);
        class{end+1} = i;
        [r,~]=find(ismember(CP,i));
        while length(unique(reshape(CP(r,:),[],1))) > length(class{end})
            class{end} = unique(reshape(CP(r,:),[],1));
            unclassified = setdiff(unclassified,class{end});
            i = class{end};
            [r,~]=find(ismember(CP,i));
        end
    end
    

    % *Group mutually concordant complexes with balanced*
   clear CP
   [CP(:,1),CP(:,2)] = find(CC~=0);
   unclassified = unique(CP);
    
    class_with_balanced=[];
    while ~isempty(unclassified)
        i = unclassified(1);
        class_with_balanced{end+1} = i;
        [r,~]=find(ismember(CP,i));
        while length(unique(reshape(CP(r,:),[],1))) > length(class_with_balanced{end})
            class_with_balanced{end} = unique(reshape(CP(r,:),[],1));
            unclassified = setdiff(unclassified,class_with_balanced{end});
            i = class_with_balanced{end};
            [r,~]=find(ismember(CP,i));
        end
    end
         
    save(strcat('Results_concordant_original_irreversibility_considered/',name,'.mat'))
end