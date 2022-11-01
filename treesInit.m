% trees implementation
clear all;
d = uigetdir();
cd(d);
cn2 = dir(d);

cn = length(cn2)-2;

folders = [{cn2(3:end).name}']
cd ../layer-analysis
%%


addpath(['/Users/benjaminsylvanus/Documents/' ...
        'CoopMRI/treestoolbox-master/']);

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

clear t
for i = 1:length(folders)
    try 
            t{i}=load_tree(fullfile(d,folders{i}));
    catch ME
        disp(ME)
        disp("L")
    end
end




