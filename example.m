clear
addpath('functions/')
% example network
%        v1 v2 v3 v4 v5 v6 v7 v8 v9 v10
model.S=[-2  2  0  0  0  0  1 -1  0  0; % A
          1 -1 -1 -1  1  2 -2  0 -2  2; % B
          0  0  1 -1  1  0  0  0  0  0; % C
          0  0  0  1 -1 -1  0  0  0  0; % D
          0  0  0  0  0  0  1 -1  0  0; % E
          0  0  0  0  0  0  0  1  1 -1];% F
model.mets={'A';'B';'C';'D';'E';'F'};
model.rxns={'2A->B';'B->2A';'B->C';'B+C->D';'D->B+C';...
    'D->2B';'2B->A+E';'A+E->F';'2B->F';'F->2B'};
model.c=zeros(size(model.rxns));
model.lb=ones(size(model.rxns))*1e-3;
model.ub=ones(size(model.rxns))*1000;
model.b=zeros(size(model.mets));
model.csense=repmat('E',size(model.b));

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

%% species degree

species_degree_complexes = sum(model.Y~=0,2);
species_degree_reactions = sum(model.S~=0,2);

%% balanced complexes
[B,group] = find_balanced_complexes(model);
model.complexes(B) 

% .complexes includes complex names in the format 
% stoichiometry*index of metabolite
% i.e. 
% >> model.complexes
% 
% ans =
% 
%   8×1 cell array
% 
%     {'2*1'    }  --> 2A --> balanced
%     {'1*2'    }  --> B
%     {'1*3'    }  --> C
%     {'1*2+1*3'}  --> B+C
%     {'1*4'    }  --> D --> balanced
%     {'2*2'    }  --> 2B
%     {'1*1+1*5'}  --> A+E --> balanced
%     {'1*6'    }  --> F --> balanced
%      
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

