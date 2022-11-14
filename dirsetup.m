% CELLTYPE / {L1:L6} / BYID / { swcfile, statsfile }
% 
% append dirs?
% 
% CELLTYPE / {data, stats, swc} / {L1:L6} / CELLID.{ske, skv, skr}
% 
% now: 
% 
% CELLTYPE / swc / {L1:L6} / CELLID.swc;
% 
% Read from data/L{n} -> Place into swc/L{n};
% 


%%
addpath("./EMSeg/");
addpath("./src/");
addpath(genpath("./treestoolbox-master/"));
addpath(genpath("./layer-data/"));
addpath(genpath("/autofs/space/symphony_002/users/BenSylvanus/Cells"));
%% Generate Folders
clear all;
root = "/autofs/space/symphony_002/users/BenSylvanus";
CELL_FOLDER_PATH = "/Cells";

% Cell type 
%   - Data
%   - SWC
%   - STATS

cd (root+CELL_FOLDER_PATH);
d = dir();

folders = d(3:end)';

for i = 1:length(folders)
   cd(append(root,CELL_FOLDER_PATH,'/',folders(i).name));
   layer_dir = dir();
   layer_dir = layer_dir(3:end)';
   s1 = sprintf('\n\n%s', ...
       append(root,CELL_FOLDER_PATH,'/',folders(i).name));
   disp(s1);
   layers = dir('./data');

   for j = 3:length(layers)
       str = sprintf('|%s%s', "    ", ...
            append('../',"./swc/"+layers(j).name));
       disp(str);
   end
end

% SWC_FOLDER_PATH = CELL_FOLDER_PATH + TYPE 


PACKAGE_PATH = "/layer-analysis-main";
% 
% % d = uigetdir();
% 
% cd(append(root,CELL_FOLDER_PATH));
% cell_dir = dir();
% 
% folders = cell_dir(3:end)';
% 
% for i = 1:length(folders)
%    cd(append(root,CELL_FOLDER_PATH,'/',folders(i).name));
%    layer_dir = dir();
%    layer_dir = layer_dir(3:end)';
%    s1 = sprintf('\n\n%s', ...
%        append(root,CELL_FOLDER_PATH,'/',folders(i).name));
%    disp(s1);
%    for j = 1:length(layer_dir)
%         d = dir(layer_dir(j).name);
%         str = sprintf('|%s%s\tCellcount = %d', "    ", ...
%             append('../',layer_dir(j).name), floor((length(d)-2)/3));
%         disp(str);
%    end
% end

cd(root+PACKAGE_PATH);


%%
% d = uigetdir();

cd(append(root,CELL_FOLDER_PATH));
cell_dir = dir();

folders = cell_dir(3:end)';
clear s;
s = struct();

for i = 1:length(folders)
   s(i).type = folders(i).name;
   cd(append(root,CELL_FOLDER_PATH,'/',folders(i).name));
   layer_dir = dir();
   layer_dir = layer_dir(3:end)';
   for j = 1
        d = dir(layer_dir(j).name);
        cd(layer_dir(j).name);
        for k = 3:length(d)
            s(i).data(k-2).Layer = string(d(k).name);
            cd(d(k).name);
            skr = dir("skr*.txt");
            ske = dir("ske*.txt");
            skv = dir("skv*.txt");
            s(i).data(k-2).ske = string({ske.name}');
            s(i).data(k-2).skv = string({skv.name}');
            s(i).data(k-2).skr = string({skr.name}');
            if numel(skr)>0
            skv = string({skv.name}');
            ske = string({ske.name}');
            skr = string({skr.name}');
            skv = split(skv,[".txt","skv"]);
            ske = split(ske,[".txt","ske"]);
            skr = split(skr,[".txt","skr"]);
            
%             s(i).data(k-2).CID = (skv(:,2));
            
            s(i).data(k-2).CID = loadCells(ske(:,2),skv(:,2),skr(:,2));
            end
            cd ../;
        end
        cd ../;
   end
end

cd(root+PACKAGE_PATH);
%%
parfor i = 1:length(s)
    disp(s(i).type);
    ctype = s(i).type;
    for j = 1:length(s(i).data)
        layer = s(i).data(j).Layer;
        disp(s(i).data(j).Layer);
        s(i).data(j).g = {};
           for k = 1:length(s(i).data(j).CID)
               try
               cid = s(i).data(j).CID(k);
               p = append(root,CELL_FOLDER_PATH,"/",ctype,"/data/",layer);
               pswc = append(root,CELL_FOLDER_PATH,"/",ctype,"/swc/",layer);
               G = loadsk(p,cid);
               swc = fixcell_Soma(G,cid,pswc);
               s(i).data(j).g = [s(i).data(j).g; {G}];
               catch ME
                   disp(ME);
               end
           end
    end
end


%%
% for i = 1:length(layer)
%     ske = layer(i).ske;
%     skv = layer(i).skv;
%     skr = layer(i).skr;
% 
%     s = loadCells(ske,skv,skr);
% end


%% Functions

function C = loadCells(ske,skv,skr)
    C = intersect(skr,intersect(ske,skv));
end

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


