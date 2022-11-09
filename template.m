%{ 
    {root: type: layer: files:}
%}

for i = 1:length(layer)
    ske = layer(i).ske; 
    skv = layer(i).skv; 
    skr = layer(i).skr;
    s = loadCells(ske,skv,skr);
end

%% FUNCTIONS 
function s = loadCells(edge,vert,radius);
    C = intersect(radius,intersect(edge,vert));
    for i = 1:length(C)
        G = loadsk(mkpath(root,type,layer,files), C(i));
    end
end

function p = mkpath(varargin)
    p = strjoin(varargin,'/');
end

function G = loadsk(path,id)
    %{
     * Generates an Undirected Graph from Node Edge Radius files  
    Inputs:
    --***************************************************************--
    root * location of data
    id   * Cell Id   
    Outputs:
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
    s = ske(:,1) + 1; t = ske(:,2)+1;
    x = skv(:,1); y = skv(:,2); z = skv(:,3);

    Name = string(1:size(x,1))';

    NodeTable = table(x,y,z,Name,skr,'VariableNames',...
        {'x' 'y' 'z' 'Name' 'Radius'});
    
    EdgeTable = table([s t], 'VariableNames',{'EndNodes'});
    G = graph(EdgeTable,NodeTable);
end





