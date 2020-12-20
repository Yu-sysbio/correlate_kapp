minclr = [255,255,178]/255;
maxclr = [242,94,13]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];

%% Growth conditions
[kapp,b,~] = xlsread('pnas.1514240113.sd01.xlsx','kapp 1s');
rxn = b(4:end,1);
conditions = b(3,3:end);

nanidx = any(isnan(kapp),2);
rxn = rxn(~nanidx);
kapp = kapp(~nanidx,:);

clear b nanidx;

label = num2cell(1:1:length(conditions));
data = kapp;

[RHO,~] = corr(log10(data),'type','Pearson');

figure();
h = heatmap(label,label,RHO,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','k');
title(['Correlations between kapps (N = ' num2str(length(rxn)) ') under ' num2str(length(conditions)) ' conditions']);
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
% h.Colormap = copper;

set(gcf,'position',[20 20 300 300]);
set(gca,'position',[0.15 0.15 0.65 0.65]);

%% KO ALE strains
[a,b,~] = xlsread('pnas.2001562117.sd01.xlsx','Dataset_S1A_protein_ab');
kapp_list = a(:,end);
rxnid_list = b(2:end,1);
rxnid_list = strrep(rxnid_list,'_b','_reverse');
strain_list = b(2:end,3);

nanidx = isnan(kapp_list);
kapp_list = kapp_list(~nanidx);
rxnid_list = rxnid_list(~nanidx);
strain_list = strain_list(~nanidx);

rxnlist = cell(0,1);
strainlist = unique(strain_list)';
kapplist = zeros(0,length(strainlist));

unqrxns = unique(rxnid_list);

for i = 1:length(unqrxns)
    rxn_tmp = unqrxns(i);
    idx_tmp = ismember(rxnid_list,rxn_tmp);
    strain_tmp = strain_list(idx_tmp);
    if length(unique(strain_tmp)) == length(strainlist)
        kapp_tmp = kapp_list(idx_tmp);
        kapplist_tmp = zeros(1,length(strainlist));
        for j = 1:length(strainlist)
            kapplist_tmp(1,j) = mean(kapp_tmp(ismember(strain_tmp,strainlist(j))));
        end
        rxnlist = [rxnlist;rxn_tmp];
        kapplist = [kapplist;kapplist_tmp];
    end
end
clear a b kapp_list rxnid_list strain_list nanidx unqrxns i j idx_tmp kapp_tmp kapplist_tmp rxn_tmp strain_tmp;

% label = num2cell(length(conditions)+1:1:length(conditions)+length(strainlist));
label = strainlist;
data = kapplist;

[RHO,~] = corr(log10(data),'type','Pearson');

figure();
h = heatmap(label,label,RHO,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','k');
title(['Correlations between kapps (N = ' num2str(length(rxnlist)) ') of ' num2str(length(strainlist)) ' strains']);
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';

set(gcf,'position',[320 320 220 220]);
set(gca,'position',[0.15 0.15 0.6 0.6]);

%% Combined dataset
% rxn_combined = intersect(rxn,rxnlist);
% [~,m] = ismember(rxn_combined,rxn);
% [~,n] = ismember(rxn_combined,rxnlist);
% kapp_combined = [kapp(m,:) kapplist(n,:)];
% label_combined = [conditions strainlist];
% 
% label = num2cell(1:1:length(label_combined));
% data = kapp_combined;
% 
% [RHO,~] = corr(log10(data),'type','Pearson');
% % [RHO,~] = corr(data,'type','Pearson');
% 
% figure();
% h = heatmap(label,label,RHO,'Colormap',clrmap,...
%     'ColorMethod','count','CellLabelColor','k');
% title(['Correlations between kapps (N = ' num2str(length(rxn_combined)) ') of combined datasets']);
% set(h,'FontSize',6,'FontName','Helvetica');
% h.FontColor = 'k';


