%{

folder setup:

CELLTYPE / {L1:L6} / CELLID.{ske, skv, skr}

create directory setup

store swc in 

CELLTYPE / {L1:L6} / BYID / { swcfile, statsfile }

append dirs?

CELLTYPE / {data, stats, swc} / {L1:L6} / CELLID.{ske, skv, skr}

now: 

CELLTYPE / swc / {L1:L6} / CELLID.swc;

Read from data/L{n} -> Place into swc/L{n};

Cell: 
key: id;
type: celltype;
swcfile: Cell{Key}.swc;
stats: matfile;
%}
%%
addpath("./EMSeg/");
addpath("./src/");
addpath(genpath("./treestoolbox-master/"));
addpath(genpath("./layer-data/"));
addpath(genpath("/autofs/space/symphony_002/users/BenSylvanus/Cells"));
%% Folders
clear all;

root = "/autofs/space/symphony_002/users/BenSylvanus";
CELL_FOLDER_PATH = "/Cells";
PACKAGE_PATH = "/layer-analysis-main";


% Cell type 
%   - Data
%   - SWC
%   - STATS

cd (root+CELL_FOLDER_PATH);
d = dir();
folders = d(3:end)';
%%
clear s;
s = struct();

for i = 1:length(folders)
   s(i).type = folders(i).name;
   cd(append(root,CELL_FOLDER_PATH,'/',folders(i).name));
   layer_dir = dir();
   layer_dir = layer_dir(3:end)';
   s1 = sprintf('\n\n%s', ...
       append(root,CELL_FOLDER_PATH,'/',folders(i).name));
   disp(s1);
   swc_path = dir('./swc');
   cd("./swc")
   for j = 3:length(swc_path)
       s(i).layer(j-2).Name = swc_path(j).name;
       swc_dir = dir(swc_path(j).name);
       [ps, fs, full] = extractInfo(swc_dir);
       s(i).layer(j-2).files = [ps,fs];
       s(i).layer(j-2).fullpath = full;
       str = sprintf('%s', "    ", ...
            append("./swc/"+swc_path(j).name));
       disp(str);
   end
end

% SWC_FOLDER_PATH = CELL_FOLDER_PATH + TYPE 

cd(root+PACKAGE_PATH);
%%
start_trees;
%%
% for i = 2:length(s)
for i = 1:length(s)
    disp(s(i).type);
    ctype = s(i).type;
    trees{i} = struct();
    trees{i}.type = ctype;
    trees{i}.layers = struct();
    t_i = trees{i};
    for j = 1:length(s(i).layer)
        layer = s(i).layer(j);
        trees{i}.layers(j).name = layer.Name;
        trees{i}.layers(j).t = {};
    end
end
%%
for i = 1:length(s)
    tic;
    disp(s(i).type);
    ctype = s(i).type;
    t_i = trees{i};
    l_j = t_i.layers;
    parfor j = 1:length(s(i).layer)
        tic;
        fullpath = s(i).layer(j).fullpath;
        tree_temp = {};
        for k = 1:size(fullpath,1)
            tic;
            [t, name, path] = load_tree (fullpath{k}, ' ');
            tree_temp{k} = t;
            t2 = toc;
            fprintf('\t\t-----%d:%d \t %d\n', k,size(fullpath,1),t2);
        end
        time = toc;
        str = sprintf('Completed %s/%s \t %4dseconds\n', ctype,s(i).layer(j).Name,time);
        disp(str);
        l_j(j).t = tree_temp;
    end
    t_i.layers = l_j;
    trees{i} = t_i;
end
% [tree, name, path] = load_tree (name, options)
% ----------------------------------------------
%
% Loads the metrics and the corresponding directed adjacency matrix to
% create a tree in the trees structure.
%
% Input
% -----
% - name     ::string: name of the file to be loaded, incl. the extension.
%     {DEFAULT: open gui fileselect, replaces format entry}
%     formats are file extensions:
%     '.mtr' : TREES toolbox internal format
%        (this is just a matlab workspace!)
%        such a file can contain more than one tree, up to 2 depth for e.g.
%        cgui_tree: {{treei1, treei2,... }, {treej1, treej2,...}, ...}
%        or: {tree1, tree2, ...} or just tree.
%     '.swc' : from the Neurolucida swc format ([inode R X Y Z D/2 idpar])
%        comments prefixed with "#", otherwise only pure ASCII data
%        May contain multiple trees!!
%     '.neu' : from NEURON transfer format .neu (see neu_tree)
%        not every NEURON hoc-file represents a correct graph, read
%        about restrictions in the documentation.
%     {DEFAULT: '.mtr'}
% - options  ::string:
%     '-s'   : show
%     '-r'   : repair tree, preparing trees for most TREES functions
%     '-ks'  : keep sections from NEURON as regions
%     {DEFAULT: '-r' for .swc/.neu 'none' for .mtr}
%
% Output
% ------
% If no output is declared the tree is added to trees
% - tree     :: structured output tree
% - name     ::string: name of output file;
%     []     no file was selected -> no output
% - path     ::string: path of the file
%   complete file name is therefore: [path name]


%% Stats Loop

astrocyte = trees{1}.layers;
interneuron = trees{2}.layers;
microglia_opc = trees{3}.layers;
oligodendrocyte = trees{4}.layers;
pyramidal = trees{5}.layers;
spiny_stellate = trees{6}.layers;
%%

%%
stats = {}; 

parfor i = 1:length(trees)
    stats{i} = stats_tree({trees{i}.layers.t},{trees{i}.layers.name},[],'-x');
end

%%
tic;
clear spiny_stats;
tree = {spiny_stellate(:).t};
name = {spiny_stellate(:).name};
spiny_stats = stats_tree(tree,name,[],'-x');
toc;

tic;
clear pyramidal_stats;
tree = {pyramidal(:).t};
name = {pyramidal(:).name};
pyramidal_stats = stats_tree(tree,name,[],'-x');
toc;

tic;
clear astro_stats;
tree = {astrocyte(:).t};
name = {astrocyte(:).name};
astro_stats = stats_tree(tree,name,[],'-x');
toc
tic;
clear interneuron_stats;
tree = {interneuron(:).t};
name = {interneuron(:).name};
interneuron_stats = stats_tree(tree,name,[],'-x');
toc;

clear microglia_stats;
tree = {microglia_opc(1:length(microglia_opc)-1).t};
name = {microglia_opc(1:length(microglia_opc)-1).name};
microglia_stats = stats_tree(tree,name,[],'-x');

tic;
clear oligodendrocyte_stats;
tree = {oligodendrocyte(:).t};
name = {oligodendrocyte(:).name};
oligodendrocyte_stats = stats_tree(tree,name,[],'-x');
toc;
%%
close all;

figure(Name="Astrocyte Stats");
dstats_tree(astro_stats);

figure(Name="Interneuron Stats");
dstats_tree(interneuron_stats);

figure(Name="Spiny Stats");
dstats_tree(spiny_stats);

figure(Name="Pyramidal Stats");
dstats_tree(pyramidal_stats)

figure(Name="Microglia Stats");
dstats_tree(microglia_stats);

figure(Name="Oligodendrocyte Stats");
dstats_tree(oligodendrocyte_stats);


%%
for i = 1:length(stats)
    figure(i);
    dstats_tree(stats{i});
end
%% Histogram

comb_stats = [astro_stats, interneuron_stats, spiny_stats, ...
    pyramidal_stats, microglia_stats, oligodendrocyte_stats];
comb_title = ["astro_stats" "interneuron_stats" "spiny_stats" ...
    "pyramidal_stats" "microglia_stats" "oligodendrocyte_stats"];

%%

%{
vol fraction w/o std;
bar plot for each layer by cell type (include std);
%}




%%
%%
close all;
hold on;
for j = 1:length(comb_stats)
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

%%
[h,q] = plotbylayer(stats,order);
%% Functions
function s = extractStats(l,stats,order)
    s = struct();
    s.dstats = struct();
    s.gstats = struct();
    s.s = {};
    for i = 1:length(stats)
        idx = stats{i}.s == l;
%         s(i).type = order(i);
        s.dstats(i).BO = stats{i}.dstats(idx).BO;
        s.dstats(i).Plen = stats{i}.dstats(idx).Plen;
        s.dstats(i).peucl = stats{i}.dstats(idx).peucl;
        s.dstats(i).angleB = stats{i}.dstats(idx).angleB;
        s.dstats(i).blen = stats{i}.dstats(idx).blen;



%         s.gstats(i) = stats{i}.gstats(idx);
%         s.s{i} = order(i);
%         s.dim = stats{i}.dim;
    end
end

function [file, path, fullpath] = extractInfo(S)
    elements = S(3:end);
    file = {elements.name}';
    path = {elements.folder}';
    fullpath = fullfile(path,file);
end

function [h,q] = plotbylayer(stats,order)
    Layers = ["L1","L2","L3","L4","L5","L6"];

    for i = 1:length(Layers)
        q(i).Layer = Layers(i);
        q(i).s = extractStats(Layers(i),stats,order);
    end
    h = 0;
end



function order = getorder(trees)
    for i = 1:length(trees)
        order(i) = string(trees{i}.type);
    end
end


