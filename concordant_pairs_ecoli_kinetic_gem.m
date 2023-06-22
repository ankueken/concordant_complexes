clear
addpath('functions/')
% e. coli kinetic network

load('../Ecoli_kinetic_models/genome_scale/Data.mat')

model.S = network_data.S_f_b;
model.b = zeros(size(model.S,1),1);
model.c = zeros(msize(model.S,2),1);
model.mets = network_data.metab;
model.rxns = network_data.rxn_f_b;
model.lb = zeros(size(model.S,2),1);
model.ub = ones(size(model.S,2),1)*1000;
model.rxnNames = network_data.rxn_f_b;
model.metNames = network_data.metab;
model.rev = zeros(size(model.lb));
model.csense = repmat('E',size(model.b));

% get A and Y matrix
save('model_temp.mat')
cd functions
system(strjoin({'"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"', 'get_AY_matrix.r','../model_temp.mat'}));
cd ..
A=importdata('model_temp.A');
model.A=A.data;
model.complexes=A.textdata(2:end);
Y=importdata('model_temp.Y');
model.Y=Y.data;
clear A Y
delete('model_temp.A','model_temp.Y','model_temp.mat')
Results_balanced.MODEL_r{1}=model;
%% species degree

species_degree_complexes = sum(model.Y~=0,2);
species_degree_reactions = sum(model.S~=0,2);

%% balanced complexes
[B,group] = find_balanced_complexes(model);
model.complexes(B) 

%% 1 - mutually concordant balanced complexes
 
CC = zeros(length(model.complexes));
CC(B,B) = 1;

%% -1 - trivial concordant complexes
% complexes including degree 2 species

idx = find(species_degree_complexes==2);
for i=1:length(idx)
    tcc = find(model.Y(idx(i),:)~=0);
    if CC(tcc,tcc)~=1
        CC(tcc,tcc) = -1;
        CC(logical(eye(size(CC)))) = 0;
    end
end

%% 2 - other couplings

coupling_pairs = find_coupling_complexes(model,B,1,0,group);

for i = 1:size(coupling_pairs,1)
    if CC(coupling_pairs(i,1),coupling_pairs(i,2)) == 0
        CC(coupling_pairs(i,1),coupling_pairs(i,2)) = 2;
        CC(coupling_pairs(i,2),coupling_pairs(i,1)) = 2;
    end
end

%% concordance modules

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
    
%% plot

CC_triu = triu(CC);

balanced_concordant = sum(sum(CC_triu==1));
trivial_concordant = sum(sum(CC_triu==-1));
higher_order_concordant_pairs = sum(sum(CC_triu==2));
total=balanced_concordant+trivial_concordant+higher_order_concordant_pairs;

bar(([balanced_concordant trivial_concordant higher_order_concordant_pairs]/total*100))
ylabel('Concordant complex pairs (%)')
set(gca,'XTickLabels',1:3,'XTickLabels',{'balanced' 'trivial' 'higher order'})

