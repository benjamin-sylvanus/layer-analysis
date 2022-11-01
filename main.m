
%%
clear all;
d = uigetdir();
cd(d);
cn2 = dir(d);

cn = length(cn2)-2;

folders = [{cn2(4:end).name}']
cd ../layer-analysis

% palombo paper figure out soma representation.


%%
addpath('/Users/benjaminsylvanus/Documents/GitHub/FixCell/New FC/lib/EMSeg');
count = [];

for i=1:numel(folders)

    try
    id = string(folders{i});
    path =d+"/"+id;
    G = loadsk(path, id);
        close all;
        p = append('/Users/benjaminsylvanus/Desktop/5-swc-sf/',id,'/');
        if isfolder(p)
            cd(p);
            delete *.swc *.mat;
            rmdir(p);
            cd ('/Users/benjaminsylvanus/Desktop/layer-analysis/')
        end
        mkdir(p);
        [swc] = fixcell_Soma(G,id,p);

        count(i,1) = length(G.Nodes.x); 
        count(i,2) = swc.mt; count(i,3) = id;
        save(p+'values.mat', 'count');
    end
end


%%
function h = plotsk2(G)
            %{
        * Plots the graph on 3D axis
        
        Inputs:
        --***************************************************************--
        G * Graph to plot
        c * Node Color
        
        OUT:
        --***************************************************************--
        h * Axis Handle
            %}
            h = plot(G,'XData',G.Nodes.x,'YData',G.Nodes.y,'ZData',G.Nodes.z);
            h.NodeLabel = {};
end
%%

%{

Generic outline

for i in folders

1. if folder i contains .swc
continue

else make swc

2. if folder contains config.txt
continue 
else 
analyze swc.

%}

%%

% function [ske, skv, skr, G, x, y, z, Name] = loadsk(root,id) 

function [G] = loadsk(path,id)
%{
 * Generates an Undirected Graph from Node Edge Radius files
    
Inputs:
--***************************************************************--
root * location of data
    
 OUT:
--***************************************************************--
G * the undirected graph
    * Unused Outputs:
    ske * edge data from file
    skv * vertex data from file
    skr * radius data from file
    x * x data 
    y * y data 
    z * z data 
    Name * node names
%}


    ske = load(fullfile(path,sprintf('ske%s.txt',id)));
    skv = load(fullfile(path,sprintf('skv%s.txt',id)));
    skr = load(fullfile(path,sprintf('skr%s.txt',id)));
   
    s = ske(:,1)+1;
    t = ske(:,2)+1;
    
    x = skv(:,1);
    y = skv(:,2);
    z = skv(:,3);
    
%     x = x-min(x);
%     y = y-min(y);
%     z = z-min(z);
    
    Name = string(1:size(x,1))';
    
    NodeTable = table(x,y,z,Name,skr,'VariableNames',...
        {'x' 'y' 'z' 'Name' 'Radius'});
    
    EdgeTable = table([s t], 'VariableNames',{'EndNodes'});
    G = graph(EdgeTable,NodeTable);
end


