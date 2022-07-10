%% Script to Compute Pathway Score (PC1) for curated KEGG Pathways 

cd /Users/sarahcherkaoui/MetabCCLs/data

% Metabolomics Data
load('metabolomicsData.mat')
% KEGG Pathways
load('modelHSA.mat')
 
%To Merge replicate - cell lines 
[uniqueCell,ib,ic] = unique(fe3.dsSamplePerturbation);
fe3.dsUser2=fe3.dsSamplePerturbation

% Log2 Zscore
fe3.merged.mat=zscore(log2(fe3.merged.mat)')'

% Mapping to Pathways
fe3.merged.Pathways=modelHSA.Pathways;
fe3.merged.PathwaysName=modelHSA.PathwaysName;
fe3.merged.IonToPathway=zeros(length(fe3.merged.Pathways),length(fe3.merged.ionTopId));

for i=1:length(fe3.merged.ionTopId)
    posMetab=find(ismember(modelHSA.Metabolites,strsplit(char(fe3.merged.ionTopId(i)),'; ')));
    for j=1:length(posMetab)
        posPath=find(modelHSA.MetaboliteToPathway(:,posMetab(j)));
        fe3.merged.IonToPathway(posPath,i)=1;
    end
end


%% 
% Generate PCs for all Pathways
%%
typesPath=cell(3,length(fe3.merged.Pathways)); 
typesPathPCAScaled=[];
typesLoadings=cell(3,length(fe3.merged.Pathways)); % save the loadings / contribution of each metabolites to the PC1
typesLoadingsScaled=cell(length(fe3.merged.Pathways),1);
explainedVar=[];
explainedVar2=[];

PathwUsed=[];

plotting = 1;
for i=1:length(fe3.merged.Pathways)
    metabolites=find(fe3.merged.IonToPathway(i,:));
    if(length(metabolites)>=4) %No Ox Phos
        dataMatrix = fe3.merged.mat(metabolites,:)';

        [coeff,score,latent,tsquared,explained,mu] = pca(dataMatrix,'NumComponents',3);
        PCAcoor = dataMatrix*coeff;
        clear title xlabel ylabel

        explainedVar=[explainedVar explained(1)];
        explainedVar2=[explainedVar2 explained(2)];

        [p,d] = size(coeff);
        [~,maxind] = max(abs(coeff),[],1);
        colsign = sign(coeff(maxind + (0:p:(d-1)*p)));
        coef = bsxfun(@times,coeff,colsign);
        maxCoefLen = sqrt(max(sum(coeff.^2,2)));
        
        typesPath{1,i}=score(:,1);
        typesPath{2,i}=score(:,2);
        typesPath{3,i}=score(:,3);
        % Need to decide which Loading you wanna show 
        typesLoadings{1,i}=coeff(:,1)';
        typesLoadings{2,i}=coeff(:,2)';
        typesLoadings{3,i}=coeff(:,3)';

        scores = bsxfun(@times, maxCoefLen.*(score ./ max(abs(score(:)))), colsign);
        typesPathPCAScaled=[typesPathPCAScaled scores(:,1)]; %%%%%%%%%%%%%%%% If you want normalise score or not
        typesLoadingsScaled{i}=coef(:,1);
        PathwUsed=[PathwUsed fe3.merged.PathwaysName(i)];
    else
        typesPathPCAScaled=[typesPathPCAScaled real(nan(length(fe3.dsSamplePerturbation),1))];
    end
end

%% 
% Merge PC1 
%%
s=size(typesPathPCAScaled);
replicatesHeatmap=nan(length(uniqueCell),s(2));
replicatesSTD=nan(length(uniqueCell),s(2));

% Merge by median all biological/technical replicates & Variance
for j = 1:s(2)
    % For all cell lines
    for i = 1:length(uniqueCell)
        if(~isnan(typesPathPCAScaled(1,j))) 
            pos=find(ic==i);
            repcluster=typesPathPCAScaled(pos,j);
            % Store variance of each Cell lines to a pathways
            replicatesSTD(i,j)=std(repcluster);
            replicatesHeatmap(i,j)= median(repcluster);
        end       
    end
end

% Plot standard deviation
posNan = find(~isnan(replicatesSTD(1,:))); % position without Nan


%%
% Remove individual with extreme standard deviation
%%

removeIndividual=find(replicatesSTD>1);

% Do not remove bad pathways or cell line but replcae to 0 in replicatesHeatmap
replicatesHeatmap(removeIndividual)=0;

%% To export for R
writetable(table(PathwUsed),"Pathways.csv",'WriteVariableNames',false)
writetable(table(cellstr(uniqueCell)),"CellLines.csv",'WriteVariableNames',false)
csvwrite("PathwayScore_180CCL.csv",replicatesHeatmap(:,posNan))

% For export of Path
CellLineMeta=table(cellstr(uniqueCell)',fe3.dsUser1(ib),fe3.dsSampleAmount(ib),fe3.dsUser4(ib),fe3.dsUser3(ib),fe3.dsPlate(ib)')
writetable(CellLineMeta,"CellLinesMeta.csv",'WriteVariableNames',false)

