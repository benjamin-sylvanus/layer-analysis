
%{
Index in position 2 exceeds array bounds.

Error in vhull_tree (line 116)
X            = [(X (ipart)); (points (:, 1))] + DD (1);

Error in stats_tree (line 261)
                [~, ~, ~, vol]                 = vhull_tree ( ...

Error in treesInit (line 42)
stats = stats_tree(tnot,[],[],'-s');
%}

% trees implementation
clear all;
d = uigetdir();
cd(d);
cn2 = dir(d);
cn = length(cn2)-2;
folders = [{cn2(3:end).name}'];
cd ../

%%
addpath("./EMSeg/");
addpath("./src/");
addpath(genpath("./treestoolbox-master/"));
addpath(genpath("./layer-data/"));
%%

[filename, pathname] = uigetfile('*.swc', 'Pick a SWC code file');
if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
end
start_trees;
t = load_tree(fullfile(pathname,filename));

plot_tree(t);

stats_tree(t);

stats_tree(t);

%%
%     '-s'  : show results
%     '-w'  : waitbar
%     '-2d' : 2d tree, concerns hulls
%     '-x'  : no extras (much much less time consuming)
%     '-f'  : save as file
start_trees;

tnot = tm;
tic;
stats = stats_tree(tnot,[],[],'-x -s');
toc;


%% Get Stats
%     '-s'  : show results
%     '-w'  : waitbar
%     '-2d' : 2d tree, concerns hulls
%     '-x'  : no extras (much much less time consuming)
%     '-f'  : save as file
clear all;
tmp1 = load('t.mat');
tmp2 = load('t300_800.mat');


%%
q = tmp1.t;
p = tmp2.t300_800;

lq = length(q);
lp = length(p);
q(lq+1:lq+lp) = p;
t=q;



%%
clear all;
start_trees;
tic;
tmp = load("ct.mat");
t=tmp.ct;
stats = stats_tree(t,[],[],'-x -s');
toc;
%%


start_trees;
% clear t;

% ld2 300-800

tic;
for i = 1:length(folders)
    tic;
    disp("current: "+ i);
    disp("Loading swc...");
    
    pathname = append(d,'/',folders{i});
    filename = append(folders{i},'.swc');
    file = fullfile(pathname,filename);
    try 
            t300_800{i}=load_tree(file);
    catch ME
        disp(ME);
    end
    toc;
end
save("t300_800.mat",'t300_800');
toc;

%%
gs = stats.gstats;
figure();
% figure(3);

content = fieldnames(gs);

for i = 1:length(content)
    try
        subplot(2,7,i);
        s = sprintf("%s",content{i});
        element = gs.(s);
        histogram(element);
        title(s);
    end
end

%%







