function concordant_randomized_networks(f)
% Input: number of sample to analyze
%
% The results will be saved in folder 'Results_concordant_randomization_irreversibility_considered/'
%
% before usage update path to R in line 32:
% Results_balanced.MODEL_r{1}=clean_model(name,model,'irreversibility_considered','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"',1e-6);

addpath('functions')
mkdir('Results_concordant_randomization_irreversibility_considered/')
mkdir('Results/final_models_irreversibility_considered/')

    Data=importdata(strcat('randomization/Ec_iJO1366-matrices/Ec_iJO1366.massbalance.stmatrix.',num2str(f)));
    model.S=Data.data;
    for i=1:size(model.S,1)
        model.mets{i,1}=strcat('M',num2str(i));
    end
    model.rev=importdata('randomization/Ec_iJO1366.stmatrix.rev')';
    for i=1:size(model.S,2)
        model.rxns{i,1}=strcat('R',num2str(i));
    end
    model.lb=zeros(size(model.S,2),1);
    model.lb(model.rev==1) = -1000;
    model.ub=ones(size(model.S,2),1)*1000;
    model.c=zeros(size(model.S,2),1);
    model.c(2582) = 1;
    model.b = zeros(size(model.S,1),1);
    model.csense = repmat('E',size(model.b));

     name = strcat('Ecoli_randomization_',num2str(f));
 
        %% create models (no blocked, reaction split, A and Y matrix)
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

        coupling_pairs = find_coupling_complexes(Results_balanced.MODEL_r{1},B,f,0,group);

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

        save(strcat('Results_concordant_randomization_irreversibility_considered/',name,'.mat'))
end
