%% *Script to Compute PC1 score for all Pathways *
% *FiaMiner Pipeline Location:* Y:\users\Sarah\FiaMiner-Pipeline\AZ_Normalization_Comparison_Petra.mat
% 
% *Input* : mtt (ouput from Metabotyping function in fiaMiner- uses for mapping metabolites to pathways)
% fe2 (output from pipeline)
% KEGG (network - only use to convert )
% load('KEGG.mat')

addpath('\\imsbnas.ethz.ch\Sauer1\users\Sarah\MATLAB\Functions\gramm-master')
addpath('Y:\users\Sarah\MATLAB\FIA_New')
load('A_Data.mat')
load('modelHSA.mat')
 
 
%To Merge replicate - cell lines 
[uniqueCell,ib,ic] = unique(fe2.dsSamplePerturbation);
fe2.dsUser2=fe2.dsSamplePerturbation

%% Log2 Zscore
fe2.merged.mat=zscore(log2(fe2.merged.mat)')'

fe2.merged.Pathways=modelHSA.Pathways;
fe2.merged.PathwaysName=modelHSA.PathwaysName;
fe2.merged.IonToPathway=zeros(length(fe2.merged.Pathways),length(fe2.merged.ionTopId));

for i=1:length(fe2.merged.ionTopId)
    posMetab=find(ismember(modelHSA.Metabolites,strsplit(char(fe2.merged.ionTopId(i)),'; ')));
    for j=1:length(posMetab)
        posPath=find(modelHSA.MetaboliteToPathway(:,posMetab(j)));
        fe2.merged.IonToPathway(posPath,i)=1;
    end
end


%%
% To generate PCA for all pathways
typesPath=cell(3,length(fe2.merged.Pathways)); % save the best clustering
typesPathPCAScaled=[];
typesLoadings=cell(3,length(fe2.merged.Pathways)); % save the loadings / contribution of each metabolites to the PC1
typesLoadingsScaled=cell(length(fe2.merged.Pathways),1);
explainedVar=[];
PathwUsed=[];

plotting = 1;
for i=1:length(fe2.merged.Pathways)
    metabolites=find(fe2.merged.IonToPathway(i,:));
    if(length(metabolites)>=4) %No Ox Phos
        dataMatrix = fe2.merged.mat(metabolites,:)';

        [coeff,score,latent,tsquared,explained,mu] = pca(dataMatrix,'NumComponents',3);
        PCAcoor = dataMatrix*coeff;
        clear title xlabel ylabel
%         figure
%         explainedVar=[explainedVar explained(1)];
%         b=gscatter(PCAcoor(:,1),PCAcoor(:,2));
%         xlabel(explained(1));
%         ylabel(explained(2));
       %  if(plotting)
       %     figure
       %     a=gscatter(PCAcoor(:,1),PCAcoor(:,2), uniqueCell,'','',10)
       %     xlabel(explained(1));
        %    ylabel(explained(2));
        %    title(fe2.merged.PathwaysName(i)); 
        
        %  end

%         figure
%         hFig=histogram(PCAcoor(:,1));
%         biplot(coeff,'scores',score,'varlabels',metabolitesKEGGNames{i});
%         xlabel(explained(1));
%         ylabel(explained(2));
%         title(mtt.pwLabel(i)); 

        % To scale the scores - as seen in 
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
        
         if(plotting)
            figure
            clear g

            [typesLoad,orderLoading]=sort(typesLoadings{1,i})
            metabNames=fe2.merged.ionTopName(metabolites);
            g=gramm('x',cellstr(metabNames(orderLoading)),'y',typesLoad);
            g.stat_summary('geom',{'bar'},'dodge',0);
            g.set_title(fe2.merged.PathwaysName(i),'FontSize',12);
            g.set_order_options('x',0);
            g.axe_property('XTickLabelRotation',60);
            g.draw();

            figure
            [typesLoad,orderLoading]=sort(typesLoadings{2,i})
            metabNames=fe2.merged.ionTopName(metabolites);
            g=gramm('x',cellstr(metabNames(orderLoading)),'y',typesLoad);
            g.stat_summary('geom',{'bar'},'dodge',0);
            g.set_title(fe2.merged.PathwaysName(i),'FontSize',12);
            g.set_order_options('x',0);
            g.axe_property('XTickLabelRotation',60);
            g.draw();
        end
        
        scores = bsxfun(@times, maxCoefLen.*(score ./ max(abs(score(:)))), colsign);
        typesPathPCAScaled=[typesPathPCAScaled score(:,1)]; %%%%%%%%%%%%%%%% If you want normalise score or not
        typesLoadingsScaled{i}=coef(:,1);
        PathwUsed=[PathwUsed fe2.merged.PathwaysName(i)];
    else
     %   typesPathPCA=[typesPathPCA nan(length(fe.dsSamplePerturbation),1)];
        typesPathPCAScaled=[typesPathPCAScaled real(nan(length(fe2.dsSamplePerturbation),1))];
    end
end
%% 
% Stored PCA typing
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
clusterRep=clustergram(replicatesSTD(:,posNan)','Colormap',redbluecmap(11),'RowLabels',PathwUsed,'ColumnLabels',uniqueCell)
%% 
% We visualize result using clustergram.
% 
% ! Problem with clustergram in Matlab. Doesn't accept NaN values and scales 
% in values - Use R instead
% 
% Remove individual with extreme standard deviation
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGED TO 0.5 Because we use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normalized PCA
removeIndividual=find(replicatesSTD>1);

% Do not remove bad pathways or cell line but replcae to 0 in replicatesHeatmap
replicatesHeatmap(removeIndividual)=NaN;

% Only for plotting with clusergram- put NaN to 0 
replicatesHeatmap2=replicatesHeatmap;
replicatesHeatmap2(removeIndividual)=0;

clusterRep=clustergram(replicatesHeatmap2(:,posNan)','Colormap',redbluecmap(11),'RowLabels',PathwUsed,'ColumnLabels',cellstr(uniqueCell),'Linkage','ward')
%%
% To find cell lines with differential Warburg effect - Comparing to Shlomi Paper
[~,~,posCell]=intersect(["T47D","MCF7","BT549","HS578T","K562","CAKI1","SKOV3","TK10","MDAMB321"],uniqueCell)
 clusterRep=clustergram(replicatesHeatmap2(posCell,posNan(1:3))','Colormap',redbluecmap(11),'RowLabels',mtt.pwLabel(posNan(1:3)),'ColumnLabels',cellstr(uniqueCell(posCell)),'Linkage','ward')
[repA,repB]=sort(replicatesHeatmap2(:,posNan(2)))
HeatMap(repA','Colormap',redbluecmap(11),'RowLabels',mtt.pwLabel(posNan(2)),'ColumnLabels',cellstr(uniqueCell(repB)))

%% Save replicatesHeatmap
save replicatesHeatmap.mat replicatesHeatmap

%% To export for R
writetable(table(PathwUsed),"Pathways.csv",'WriteVariableNames',false)
writetable(table(cellstr(uniqueCell)),"CellLines.csv",'WriteVariableNames',false)
csvwrite("PCAHeatMap.csv",replicatesHeatmap(:,posNan))

% For export of Path
CellLineMeta=table(cellstr(uniqueCell)',fe2.dsUser1(ib),fe2.dsSampleAmount(ib),fe2.dsUser4(ib),fe2.dsUser3(ib),fe2.dsPlate(ib)')
writetable(CellLineMeta,"CellLinesMeta.csv",'WriteVariableNames',false)


