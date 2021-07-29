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
end

module_system_matrix=(class_count./repmat(total,size(class_count,1),1))*100;

% for which complexes at least 50% of complexes in a subsystem fall in one
% module
[r,c]=find(module_system_matrix>=50);
uniqueSubSystems(unique(c))

%%
check=find(sum(module_system_matrix)>=50);
idx=[];
for i=1:length(check)
    val=sort(module_system_matrix(:,check(i)),'descend');
    idx(i,1)=1;
    while sum(val(1:idx(i,1)))<50
        idx(i,1)=idx(i,1)+1;
    end
end

% subsystems for which number of modules in second column covers at least
% 50% of all complexes in that subsystem
table(uniqueSubSystems(check),idx);

%% figure S2
Number_modules=sum(class_count>0);
[~,idx]=sort(Number_modules);

idx(Number_modules(idx)==0)=[];

hbar=bar(module_system_matrix(:,idx)','stacked','FaceColor',[.9 .9 .9])
set(gca,'XTick',[],'XTickLabels',[])
ylabel([{'Percentage of complexes'}; {'per concordance module'}])

x = 1:length(idx);
y = sum(module_system_matrix(:,idx));
ygap = 0.1;  % Specify vertical gap between the bar and label

for i = 1:length(x) % Loop over each bar 
            xpos = x(i);        % Set x position for the text label
            ypos = y(i) + ygap; % Set y position, including gap
            htext = text(xpos,ypos,num2str(Number_modules(idx(i))));          % Add text label
            set(htext,'VerticalAlignment','bottom',...  % Adjust properties
                      'HorizontalAlignment','center')
end
ylim([0 110])  

figure
% subplot(2,1,2)
bar(Number_modules(idx)./total(idx),'FaceColor',[.9 .9 .9])
set(gca,'XTick',1:length(idx),'XTickLabels',uniqueSubSystems(idx),'XTickLabelRotation',45)   
ylabel([{'Fraction number of modules to'}; {'number of complexes per subsystem'}])
