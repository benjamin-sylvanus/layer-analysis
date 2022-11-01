classdef fixcell_Soma < handle
    
    properties (GetAccess = public, SetAccess = public)
        
        G; % Unmodified Graph
        
        filename; % Input Filename
        
        show; % - Option * Not Implemented Yet
        
        H; % Modified Graph
        
        T; % Table of Modified Graph
        
        modif; % Scaler For Writer
        
        pathin; % Path to Save Directory

        mt;
        
    end
    
    
    properties (GetAccess = private, SetAccess = private)
        
        
        typeIn; % Filetype In: * Graph || SWC *
        
        typeOut; % Filetype Out: * SWC *
        
    end
    
    methods (Access = public)
        % Constructor for fixcell_Soma
        % Takes filename or graph
        function this = fixcell_Soma(G,filename,pathin)
                this.filename = filename;
                this.pathin = pathin;
                this.typeIn = 0;

                this.G = G;
            % Set the default modifier
            this.modif = 1;
            
            % Enter Main Loop
            this.main();
        end
        
        function SortTree(this,H)
            G1 = this.H;
            
            % Reordering
            center = 1;
            
            % Find the terminal nodes in dendrites
            target = string(G1.Nodes.Name(degree(G1) == 1));
            
            % Create a directed graph from center -> Terminal Nodes
            TR = shortestpathtree(G1,string(center));
            
            % Segment ids with is soma and not soma
            idnotSoma = find(TR.Nodes.Type ~= 1);
            
            idSoma = find(TR.Nodes.Type == 1);
            ids = successors(TR,1);
            b = findnode(TR,ids);
            
            % Soma has type 1 all else is type 3
            TR.Nodes.Type(idnotSoma) = 3; TR.Nodes.Type(idSoma) = 1;
%             TR.Nodes.Type(b) = 1;
            
            
            % Sort the nodes from 1:3
            [avar,order] = sort(TR.Nodes.Type,'descend');
            
            % Reorder According to order
            orderedSTBL = reordernodes(TR,order);
            
            % Topological ordering of sorted graph
            [N1,HFast] = toposort(orderedSTBL);
            
            % Reorder with toposort graph
            orderedSTBL = reordernodes(orderedSTBL,N1);
            
            G2 = orderedSTBL;
            
            % Remove Node Names from the NodeTable
            nodes = G2.Nodes(1,[1 2 3 5 6]);
            
            % Set a temp Variable to a string of 1:number_of_nodes
            G2.Nodes.Names = string(1:numel(G2.Nodes.x))';
            
            % Set Node Names = temp Variable
            G2.Nodes.Name = G2.Nodes.Names;
            
            % Finds parent nodes for all nodes
            G3 = this.Parents(G2);
            
            % Restructure Graph to include [xyz name radius type parent]
            G4 = digraph(G3.Edges,G3.Nodes(:,[1 2 3 4 5 6 8]));
            
            % Create a table for writing purposes
            T = this.maketable(G4);
            
            this.T = T;
            
            % Write an SWC file
            this.writefile(T,'ANewSomaFixed.swc');
        end
        
        function G = Parents(this,G)
            Edges = str2double(string(G.Edges.EndNodes));
            Sorted = sortrows(Edges,2);
            temp = string(Sorted);
            ParentList = [-1 1;temp];
            G.Nodes.Parent = ParentList(:,1);
        end
        
        function main(this)
            % Determine the whether to read a graph or SWC
            typeIn = this.typeIn;
            
            if typeIn == 0
                fprintf('This will execute Next steps \n');
            end
            
            if this.typeIn == 1
                fprintf('Reading from %s \n',this.filename);
                this.readfile(this.filename);
                
            end

            
            [this,c] = this.RemoveSegments();

            if c== 1
                this.Segment;
                this.SortTree;
            else
                return;
            end
        
        end

        function [this,c] = RemoveSegments(this)
            c=0;
            G = this.G;
            [bins, binsizes] = this.conncomp(G);
            idx = binsizes(bins) == max(binsizes);
            Gp = subgraph(G, idx);
            this.mt = length(Gp.Nodes.x);
            if ((length(Gp.Nodes.x)/length(G.Nodes.x)) > 0.85)
                disp('Continuing');
                this.G = Gp;
                c=1;
            else 
                disp('Too many segmented nodes: ')
                s = sprintf( ...
                    ['Main: %d, Full %d'], ...
                    length(Gp.Nodes.x),length(G.Nodes.x));
                disp(s);
                c=0;
            end
        end

        function [bins, binsizes] = conncomp(this, G)
        %{
        * Creates bin ids for each segment in G
            * Sorted by binsize
        
        Inputs:
        --***************************************************************--
        G * the full graph
            
        OUT:
        --***************************************************************--
        bins * array that designates each node to bin [x]
        binsizes * number of nodes in each bin
        %}
            [bins0, binsizes] = conncomp(G);
            [binsizes, I]= sort(binsizes,'descend');
            bins = zeros(size(bins0));
            for i = 1:numel(binsizes)
                bins(bins0==I(i)) = i;
            end
        end
        
        function readfile(this,tempfilename)
            % Add the lib to the path
            addpath(genpath('./lib'));
            
            swcfid       = fopen ([append(this.pathin,'/',this.filename)]);
            A            = textscan (swcfid, '%s', 'delimiter', '\n');
            A            = A{1};
            fclose       (swcfid);
            
            swc          = [];
            for counter  = 1 : length (A)
                if ~isempty (A{counter})  % allow empty lines in between
                    % allow comments: lines starting with #:
                    if ~strcmp (A{counter} (1), '#')
                        swc0   = textscan (A{counter}, '%f')';
                        swc    = [swc; swc0{1}'];
                    end
                end
            end
            
            %{
                Overall Goal: Convert Soma Nodes to several nodes
                - Needs file reader
                --- Create a SWC File Reader, IMPL trees code

                --- Read into matrix and elim commented sections

                --- # inode R X Y Z D/2 idpar
            %}
            T = table(swc(:,1),swc(:,2),swc(:,3),...
                swc(:,4),swc(:,5),swc(:,6),swc(:,7),...
                'VariableNames',...
                {'Identifier',...
                'Node',...
                'XPos',...
                'YPos',...
                'ZPos',...
                'Radius',...
                'Parent'});
            
            parents = T.Parent;
            
            Edges = [parents(2:numel(parents)),...
                T.Identifier(2:numel(parents))];
            
            Nodes = parents;
            
            Name = string(T.Identifier);
            
            NodeTable = table(T.XPos,T.YPos,T.ZPos,Name,T.Radius,...
                'VariableNames',...
                {'x' 'y' 'z' 'Name' 'Radius'});
            
            EdgeTable = table(Edges, 'VariableNames',{'EndNodes'});
            
            this.G = graph(EdgeTable,NodeTable);
        end
        
        function test_read(this,a_filename)
            setfiletype(this,a_filename);
        end
    end
    
    methods (Access = private)
        
       %{ 
        function segment(this)
            G = this.G
            % Assume EMSeg on path
            % Create a Gaussian Distribution of
            %       Radius to find "Soma" nodes
            
            [mask,mu,v,p] = EMSeg([G.Nodes.Radius],6);
            
            % Create a subgraph of the nodes considered Soma
            sg = subgraph(G,mask>6);
            
            % Create Logical Arrays for Subgraph:
            %       degree <= 1 || degree > 1
            idx1 = degree(sg) <= 1;
            idx2 = degree(sg) > 1;
            
            [H,Connections] = this.saveIds();
            [X Y Z] = COM(sg);
            
            Name = string('soma');
            
            NodeTable = table(X,Y,Z,Name,max(RS.sg.Nodes.Radius),...
                'VariableNames',...
                {'x' 'y' 'z' 'Name' 'Radius'});
            
%             G1 = addnode(H,NodeTable); G1 = addedge(G1,'soma',k);
            
            
        end
%}
        
        function setfiletype(this,text)
            expression = '\.';
            extension = regexp(text,expression,'split');
            extension = extension{end};
            if extension == 'swc'
                disp("SWC")
            end
            if extension == 'txt'
                disp('Text file');
            end
        end
        
        function T = maketable(this,G)
            T = table(G.Nodes.Type,[str2double(string(G.Nodes.Name))],...
                [G.Nodes.x],[G.Nodes.y],[G.Nodes.z],...
                [G.Nodes.Radius],[str2double(string(G.Nodes.Parent))],...
                ...
                'VariableNames',...
                {'Identifier',...
                'Node',...
                'XPos',...
                'YPos',...
                'ZPos',...
                'Radius',...
                'Parent'});
        end
        
        function writefile(this,T,filename)
            
            tname = filename;
            path = '';
            
            % Open the file name
            N = strsplit(this.filename,'.');
            N{1};
            fn = append(N{1},'.swc');
            

            % extract a sensible name from the filename string:
            swc              = [T.Node T.Identifier T.XPos.*10^-1 ...
                T.YPos.*10^-1 ...
                T.ZPos.*10^-1 T.Radius.*10^-1 T.Parent];
            
            swcfile          = fopen([append(this.pathin,'/',fn)],'w'); % open file
            
            fwrite           (swcfile, ...
                ['# First Tree Made with swc extension - ' N{1},...
                char(13), char(10)],'char');
            
            fwrite           (swcfile, ...
                ['# From Trees Toolbox: procedure "swc_tree"' ...
                'part of the TREES package',...
                char(13), char(10)], 'char');
            
            fwrite           (swcfile, ...
                ['# in MATLAB',...
                char(13), char(10)], 'char');
            
            fwrite           (swcfile, ...
                ['# A copywrite',...
                char(13), char(10)], 'char');
            
            fwrite           (swcfile, ...
                ['#',...
                char(13), char(10)], 'char');
            
            fwrite           (swcfile, ...
                ['# inode R X Y Z D/2 idpar',...
                char(13), char(10)], 'char');
            
            fprintf          (swcfile, ...
                '%d %d %12.10f %12.10f %12.10f %12.10f %d\n', swc');
            %                 '%d %d %12.8f %12.8f %12.8f %12.8f %d\n', swc');
            
            fclose           (swcfile);
        end
        
        function Segment(this)
            G = this.G;
            % Assume EMSeg on path
            % Create a Gaussian Distribution of
            %       Radius to find "Soma" nodes
            f = figure();
            [mask,mu,v,p] = EMSeg([G.Nodes.Radius],4);
            close(f);
            % Create a subgraph of the nodes considered Soma
            sg = subgraph(G,mask>=3);
%             f = figure();
%             this.plotsk2(sg);



            
%             % Create Logical Arrays for Subgraph:
%             %       degree <= 1 || degree > 1
%             idx1 = degree(sg) <= 1;
%             idx2 = degree(sg) > 1;
%             
%             Edges = string(G.Edges.EndNodes);
%             
%             % from soma
%             Base = Edges(:,1);
%             % from target
%             Target = Edges(:,2);
%             
%             sgNames = string(G.Nodes.Name(mask==3));
%             sg2Names = string(G.Nodes.Name(mask~=3));
%             
%             fromSoma = ismember(Base,sgNames);
%             toNotSoma = ismember(Target,sg2Names);
%             
%             toSoma = ismember(Target,sgNames);
%             fromNotSoma = ismember(Base,sg2Names);
%             
%             idcombined = fromSoma & toNotSoma;
%             
%             idcomb2 = toSoma & fromNotSoma;
%             
%             keep = string(G.Edges.EndNodes(idcombined,2));
%             
%             keep2 = string(G.Edges.EndNodes(idcomb2,1));
%             
%             k = unique([keep;keep2]);
%             
%             H = rmnode(G,sgNames);
%             
%             [X Y Z] = this.COM(sg);
%             
%             Name = string('soma');
%             
%             NodeTable = table(X,Y,Z,Name,2*max(sg.Nodes.Radius),...
%                 'VariableNames',...
%                 {'x' 'y' 'z' 'Name' 'Radius'});
%             
%             G1 = addnode(H,NodeTable); G1 = addedge(G1,'soma',k);
            [bins, binsizes] = this.conncomp(sg);
            sg1 = subgraph(sg,bins==1);

            sNodes = ismember(G.Nodes.Name,sg1.Nodes.Name);
            other = ~ismember(G.Nodes.Name,sg1.Nodes.Name);
            G2 = G;
            G2.Nodes.Type=ones(length(G2.Nodes.x),1);
            G2.Nodes.Type(other) = 3;
            G2.Nodes.Type(sNodes) = 1;
            
%             this.H = G1;
            this.H = G2;

%             h=this.plotsk2(G2);
%             h.NodeCData = G2.Nodes.Type;

        end
        
        function [X Y Z] = COM(this,G)
            w = G.Nodes.Radius;
            x = G.Nodes.x;
            y = G.Nodes.y;
            z = G.Nodes.z;
            mass = sum(w);
            X = (sum(w.*x))/mass;
            Y = (sum(w.*y))/mass;
            Z = (sum(w.*z))/mass;
        end

        function h = plotsk2(this,G)
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
        
    end
    
end
