%%
addpath("./EMSeg/");
addpath("./src/");
addpath(genpath("./treestoolbox-master/"));
addpath(genpath("./layer-data/"));
addpath(genpath("/autofs/space/symphony_002/users/BenSylvanus/Cells"));
clear variables;
start_trees;
load('stats.mat');
load('trees.mat');
load('order.mat');


%% Thought process
stat_copy = stats;
cell_types = [...
    {'astrocyte'},{'interneuron'}, ...
    {'microglia-opc'},{'oligodendrocyte'},...
    {'pyramidal'}, {'spiny-stellate'}];


for k = 1:6
    temp = stat_copy{k};
    temp.s = cell_types;
    countr = 0;
    for z = 1:length(stats)
        i = z - countr;
        l = sprintf("L%d",k);
        idx = stats{i}.s == l;
        if (sum(idx) >0)
            temp.dstats(i) = stats{i}.dstats(idx);
            temp.gstats(i) = stats{i}.gstats(idx);
        else
            temp.gstats(i) = [];
            temp.dstats(i) = [];
            temp.s(i) = [];
            countr = countr+1;
%             names = fieldnames(temp.dstats(i));
%             for j = 1:length(names)
%                 temp.dstats(i).(names{j}) = deal({});
%             end
%             names = fieldnames(temp.gstats(i));
%             for j = 1:length(names)
%                 temp.gstats(i).(names{j}) = deal({});
%             end

            
        end
    end
    stat_copy{k} = temp;
end

%%
for i = 1:length(stat_copy)
    figure(Name=sprintf("Layer %d",i));
    tm = {stat_copy{i}.gstats.len};
    tm = cell2mat(tm(:));
    [S, L] = bounds(tm);
    edges = [S:((L-S)/5):L];
    nx = [];
    for j = 1:length(stat_copy{i}.gstats)
        nx = [nx; histcounts(stat_copy{i}.gstats(j).len,edges)];
    end
%     [N1,~] = histcounts(stat_copy{i}.gstats(1).len,edges);
%     [N2,~] = histcounts(stat_copy{i}.gstats(2).len,edges);
%     [N3,~] = histcounts(stat_copy{i}.gstats(3).len,edges);
%     [N4,~] = histcounts(stat_copy{i}.gstats(4).len,edges);
    
    bar(nx');
    xlabel("Len");
    ylabel("Number Cells");
    legend(stat_copy{i}.s);
    xticklabels(string(floor(edges(1:end-1))) + "-" +string(ceil(edges(2:end))))
    title(sprintf("Layer: %d",i));
end

%%
close all;
% for i = 2:length(stat_copy)
for i = 2
%     figure(Name=sprintf("Layer %d",i));
    tm = stat_copy{i}.gstats;

    fn = fieldnames(tm);
    n = struct();
    for k = 1:8
        nl = [];
        nx = [];
        for z = 1:length(tm)
            nxt = tm(z).(fn{k});
            nlt = repmat({stat_copy{i}.s{z}},length(nxt),1);
            nx = [nx;nxt];
            nl = [nl;nlt]; 
        end
        n(k).type = fn{k};
        n(k).data = nx;
        n(k).labels = nl;

        

        
%         [~,~,st] = anova1(nx,nl);
% 
%         [c, ~,~,gnames] = multcompare(st);
%         title(Name=sprintf("Layer %d: %s",i,string(fn{k})));
%         f = gcf;
%         fname = sprintf("./figures/Anova1/Layer%d:%s.fig",i,string(fn{k}));
%         savefig(f,fname);
%         close all;

        
    end
%     nx = [];
%     x1 = t(1).len;
%     x2 = t(2).len;
%     x3 = t(3).len;
%     x4 = t(4).len;
% 

%     xlabel("Len");
%     ylabel("Number Cells");
%     legend(stat_copy{i}.s);
%     xticklabels(string(floor(edges(1:end-1))) + "-" +string(ceil(edges(2:end))))
%     title(sprintf("Layer: %d",i));
end


%%



%%
plotbylayer(stat_copy);
%%
for j = 1:length(stat_copy)

    gs = comb_stats(j).gstats;
    content = fieldnames(gs);
    figure(Name=comb_title(j))
    for i = 1:length(content)
        try
            subplot(7,2,i)
            hold on;
            s = sprintf("%s",content{i});
            element = gs.(s);
            histogram(element);
            title(s);
        catch ME
            clc;
            disp(ME);
        end
    end
end

%% Functions
function [h] = plotbylayer(stat_copy)
    for i = 1:length(stat_copy)
        figure(Name = sprintf("L%d",i))
        dstats_tree(stat_copy{i});
    end
end



function order = getorder(trees)
    for i = 1:length(trees)
        order(i) = string(trees{i}.type);
    end
end






