clear
load('Results_concordant_original_irreversibility_considered\YeastGEM.mat')

model=Results_balanced.MODEL_r{1};
model.subSystems = cellfun(@(x) strrep(x,',',' '),model.subSystems,'UniformOutput',false);

%% remove sce number form subsystem name
for i=1:length(model.subSystems)
    for j=1:length(model.subSystems{i})
        idx=cell2mat(regexp(model.subSystems{i}(j),'sce'));
        if ~isempty(idx)
            model.subSystems{i}{j}=model.subSystems{i}{j}(idx+8:end);
        end
    end
end
for i=1:length(model.subSystems)
    for j=1:length(model.subSystems{i})
        idx=cell2mat(regexp(model.subSystems{i}(j),'^ '));
        while ~isempty(idx)
            model.subSystems{i}{j}=model.subSystems{i}{j}(idx+1:end);
            idx=cell2mat(regexp(model.subSystems{i}(j),'^ '));
        end
    end
end

CLASS=zeros(size(model.Y,2),length(class));

for j=1:length(class)
    CLASS(class{j},j)=1;
    for i=1:length(class{j})
        rxns=find(model.A(class{j}(i),:)~=0);
        
        S={strjoin(model.subSystems{rxns(1)},'/')};
        if length(rxns)>1
            for z=2:length(rxns)
                if ~contains(S{1},strjoin(model.subSystems{rxns(z)},'/'))
                    S=strcat(S{1},'/',model.subSystems{rxns(z)});
                end
            end
        end
        ComplexSubSystems{j}(i,1)={strjoin(S,'/')};
    end
end

for j=1:size(model.Y,2)
    rxns=find(model.A(j,:)~=0);
        
        S={strjoin(model.subSystems{rxns(1)},'/')};
        if length(rxns)>1
            for z=2:length(rxns)
                if ~contains(S{1},model.subSystems{rxns(z)})
                    S=strcat(S{1},'/',model.subSystems{rxns(z)});
                end
            end
        end
        ListComplexSubSystems(j,1)={strjoin(S,'/')};
end

uniqueSubSystems = unique(ListComplexSubSystems);
uniqueSubSystems = cellfun(@(x) strsplit(x,'/'),uniqueSubSystems,'UniformOutput',false);
for i=1:length(uniqueSubSystems)
    if length(uniqueSubSystems{i})>1
        uniqueSubSystems(end+1:end+(length(uniqueSubSystems{i})-1))= uniqueSubSystems{i}(2:end);
        uniqueSubSystems(i) = uniqueSubSystems{i}(1);
    end
end
uniqueSubSystems = unique([uniqueSubSystems{:}])';

for j=1:length(uniqueSubSystems)
    total(j)=sum(contains(ListComplexSubSystems,uniqueSubSystems(j)));
    for i=1:size(CLASS,2)
        class_count(i,j)=sum(contains(ListComplexSubSystems(CLASS(:,i)~=0),uniqueSubSystems(j)));
    end
    class_count(i+1,j) = total(j) - sum(class_count(:,j)); % class of complexes not in concordance relation
end

module_system_matrix=(class_count./repmat(total,size(class_count,1),1))*100;
Number_modules=sum(module_system_matrix>0);

% Suppose you have a lot of 100 floppy disks and you know that 20 of them
% are defective. What is the probability of drawing 0 through 5 defective
% floppy disks if you select 10 at random?
% p = hygepdf(0:5,100,20,10)

% for each cell we claculate a p-value
% add class not concordant 
for i=1:size(class_count,1)
    for j=1:size(class_count,2)
        if class_count(i,j) == 0
            p_hypergeom(i,j) = nan;
        else
            p_hypergeom(i,j) = hygepdf(class_count(i,j),sum(class_count,'all'),sum(class_count(:,j)),sum(class_count(i,:)));
        end
% total balls - all reactions
% marked: sum of columns - all reactions in subsystem
% drawing sum of rows - reactions in that module
    end
end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_hypergeom,0.05,'pdep','yes');

disp('percentage significant module-pathway enrichment:')
disp('significant enrichment in at least 50% of the concordance modules that')
disp('contain complexes of that subsystem')
sum((sum(adj_p<=0.05,'omitnan')./sum(~isnan(adj_p)))>=0.5)/length(uniqueSubSystems)
disp('significant enrichment in all concordance modules that')
disp('contain complexes of that subsystem')
sum((sum(adj_p<=0.05,'omitnan')./sum(~isnan(adj_p)))==1)/length(uniqueSubSystems)
