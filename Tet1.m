classdef Tet1 < matlab.mixin.Copyable
    %Tet1 P1 Tetrahedral structured mesh rectangle
    %o = Tet1(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
    
    properties
        Connectivity
        Points
        Faces
        XC
        YC
        ZC
        edges
        Element
        ElementType
        neighs
        
        nnod
        nele
        
    end
    
    
    properties (Access = private)
        x0
        x1
        y0
        y1
        z0
        z1
        nxe
        nye
        nze

    end
    
    methods
        function o = Tet1(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
            
            o.nxe = nxe; o.nye = nye; o.nze = nze;
            o.x0=x0; o.x1=x1;  o.y0=y0; o.y1=y1;  o.z0=z0; o.z1=z1;
            nx = nxe; ny = nye; nz = nze;
            
            %% Start by creating fishbone mesh in the XY-plane
            T = fishbone1(x0,x1,y0,y1,nxe,nye);
            t = T.Connectivity;
            p = T.Points;
            
            %% Extrude upwards
            dl=(z1-z0)/nze;
            z = z0:dl:z1;
            nnod = size(p,1);
            znod = reshape( ones(nnod,1)*z, [], 1);
            xynod = repmat(p,nz+1,1);
            p = [xynod,znod];
            mt = max(t(:));
            
            nzele = size(t,1);
            prism = zeros(nzele*nz,6); % preallocate space for prism elements
            l = (0:nz)*nzele;
            c = 1;
            for iz = 0:nz-1
                lo = l(c)+1;
                up = l(c+1);
                tt = [t(:,:)+mt*iz,t(:,:)+mt*(iz+1)];
                prism(lo:up,:) = tt;
                c=c+1;
            end

            prism = unique(prism,'rows');
            
            %% Prism to tet
            tet = [prism(:,[1,2,3,5]);
                   prism(:,[1,5,6,4]); 
                   prism(:,[1,5,6,3])]; 
            
            %% tet2tri
            o.nele = size(tet,1);
            o.nnod = size(p,1);
            TR = triangulation(tet,p);
            tetraTri = zeros(4*size(tet,1),3);
            o.Element(o.nele).faces = [];
            o.Element(o.nele).edges = [];
            o.edges = TR.edges;
            for iel = 1:size(tet,1)
                eli=tet(iel,:);
                iface = [eli([1,2,3]);
                         eli([1,3,4]);
                         eli([1,2,4]);
                         eli([2,3,4])];
                tetraTri(4*iel-3:4*iel,:) = iface;
                
                o.Element(iel).faces = iface;
                o.Element(iel).edges = [eli([1,2]);
                                        eli([2,3]);
                                        eli([3,4]);
                                        eli([4,1]);];

            end

            %% Output to Tet1 properties
            o.ElementType = 'Tet1Mesh';
            o.Connectivity = tet;
            o.Points = p;
            o.Faces = tetraTri;
            o.XC = p(:,1); o.YC = p(:,2); o.ZC = p(:,3);
            o.neighs = TR.neighbors;
            
            
            
            
            
            
            
            
            
            
            
        end
        
        h = vizMesh(o,varargin);
        
        [fi, fix, fiy, fiz, vol] = baseFun(o,iel,X);
        
    end
    
end


