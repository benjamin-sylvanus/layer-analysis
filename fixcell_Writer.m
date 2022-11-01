classdef fixcell_Writer < handle
    
    properties (GetAccess = public, SetAccess = public)
        G;
        filename;
        pathin;
    end
    
    properties (GetAccess = private, SetAccess = private)
        
    end
    
    methods (Access = public)
        function this = fixcell_Writer(filename,p)
            % Undirected Graph
            this.filename = filename;
            this.pathin = p;
        end
        
        function TR = write(this,G)
            this.G = G;
            
            this.G = this.IdentifyNode();

            [this.G,c] = this.RemoveSegments();
            
            if (c==1)
                this.G = this.outTree();
                
                this.G = this.fixNodeName();
                
                this.G = this.Parents(this.G);
                
    %             this.G = this.Normalize(this.G);
                
                this.write_SWC(this.filename);
            end
                TR = this.G;
            
        end
        
        function G = Normalize(this,G)
            X = G.Nodes.x; Y = G.Nodes.y; Z = G.Nodes.z;
            X = X - min(X);
            Y = Y - min(Y);
            Z = Z - min(Z);
            G.Nodes.x = X; G.Nodes.y = Y; G.Nodes.z = Z;
        end
        
                
        function  savegraph(this)
            fn = append(this.pathin,'/','Graph');
            class(fn)
            G = this.G
            save(fn,'G');
        end
        
    end
    
    methods (Access = private)

        function [G,c] = RemoveSegments(this)
            c=0
            G = this.G
            [bins, binsizes] = this.conncomp(G);
            idx = binsizes(bins) == max(binsizes);
            Gp = subgraph(G, idx);
            if ((length(Gp.Nodes.x)/length(G.Nodes.x)) > 0.85)
                disp('Continuing');
                G = Gp;
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
        
        function write_SWC(this,filename)
            T = this.maketable();
            this.writefile(T,filename);
        end
        
        function G = IdentifyNode(this)
            G = this.G;

            f = figure();
            
            [mask,mu,v,p] = EMSeg([G.Nodes.Radius],3);

            close(f);
            
            Type = mask;
            
            Type(find(mask >=2 )) = 1;
            
            Type(find(mask == 1)) = 3;
            
            G.Nodes.Type = Type;
            
        end
        
        function T = maketable(this)
            G = this.G;
            
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
            tname = this.filename + ".swc";
            path = this.pathin;
            N = split(tname,'.');

            name             = N{1};
            
            swc              = [T.Node T.Identifier T.XPos ...
                T.YPos  ... 
                T.ZPos T.Radius T.Parent];
            fi               = append(path,tname);
            addpath(path)
            swcfile          = fopen([fi],'w'); % open file
            
            fwrite           (swcfile, ...
                ['# First Tree Made with swc extension - ' name,...
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
                '%d %d %12.8f %12.8f %12.8f %12.8f %d\n', swc');
            
            fclose           (swcfile);
        end
        
        function G = outTree(this)
            
            G = this.G;
            
            wbc = centrality(G,'betweenness');

            figure();
            h = this.plotsk2(G);

            n = numnodes(G);
            
            wbc =  2*wbc./((n-2)*(n-1));

            h.NodeCData = wbc;

            f = figure();
            
            [mask,mu,v,p] = EMSeg(wbc,2);
            
            close (f);
            
            idx = find(mask == 2);
            
            H2 = subgraph(G,idx);
            
            wbc = centrality(H2,'betweenness');
            
            id = find(wbc == max(wbc));  
            if length(id)>1
                id = id(1);
            end
            
            center = string(H2.Nodes.Name(id));
            
            target = string(G.Nodes.Name(degree(G) == 1));
            
            TR = shortestpathtree(G,string(center),target);

            figure();

            h = this.plotsk2(TR)
            
            G=TR;
            
            [avar,order] = sort(TR.Nodes.Type,'ascend');
            
            orderedSTBL = reordernodes(TR,order);

            
            
            [N1,HFast] = toposort(orderedSTBL);

            orderedSTBL = reordernodes(orderedSTBL,N1);
           
            G = orderedSTBL;
        end
        
        function G = fixNodeName(this)
              G = this.G;
              G.Nodes.Names = string(1:numel(G.Nodes.x))';
              G.Nodes.Name = G.Nodes.Names;
              G = digraph(G.Edges,G.Nodes(:,[1 2 3 4 5 6 7]));

        end
        
       function G = Parents(this,G)
            Edges = str2double(string(G.Edges.EndNodes));
            Sorted = sortrows(Edges,2);
            temp = string(Sorted);
            ParentList = [-1 1;temp];
            G.Nodes.Parent = ParentList(:,1);
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



