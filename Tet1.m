classdef Tet1 < matlab.mixin.Copyable
    %Tet1 P1 Tetrahedral structured mesh rectangle
    %o = Tet1(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
    %o = Tet1(x0,x1,nxe,y0,y1,nye,z0,z1,nze,type)
    % type = 'default'|'fishbone'
    
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
        function o = Tet1(x0,x1,nxe,y0,y1,nye,z0,z1,nze,varargin)
            
            o.nxe = nxe; o.nye = nye; o.nze = nze;
            o.x0=x0; o.x1=x1;  o.y0=y0; o.y1=y1;  o.z0=z0; o.z1=z1;
            nx = nxe; ny = nye; nz = nze;
            
            %% Start by creating  mesh in the XY-plane
            if nargin == 10
                if strcmpi(varargin{1},'fishbone')
                    T = Mesh.fishbone2DP1(x0,x1,y0,y1,nxe,nye);
                    t = T.Connectivity;
                    p = T.Points;
                    
                    [elements,edgess,~]=tritoel(p(:,1),p(:,2),t);
                    tri=trinode(edgess,elements);
                    edgess=direction(edgess,p(:,1),p(:,2));
                    [tet,x3d,y3d,z3d]=extrude3d(edgess,elements,tri,p(:,1),p(:,2),z0,z1,nz);
                    p = [x3d,y3d,z3d];
                else
                    [X,Y]=meshgrid(linspace(x0,x1,nx+1),linspace(y0,y1,ny+1));
                    X = [X(:),Y(:)];
                    t = delaunay(X(:,1),X(:,2));
                    p = X;
                    
                    [elements,edgess,~]=tritoel(p(:,1),p(:,2),t);
                    tri=trinode(edgess,elements);
                    edgess=direction(edgess,p(:,1),p(:,2));
                    [tet,x3d,y3d,z3d]=extrude3d(edgess,elements,tri,p(:,1),p(:,2),z0,z1,nz);
                    p = [x3d,y3d,z3d];
                end
            else
                [X,Y]=meshgrid(linspace(x0,x1,nx+1),linspace(y0,y1,ny+1));
                X = [X(:),Y(:)];
                t = delaunay(X(:,1),X(:,2));
                p = X;
                
                [elements,edgess,~]=tritoel(p(:,1),p(:,2),t);
                tri=trinode(edgess,elements);
                edgess=direction(edgess,p(:,1),p(:,2));
                [tet,x3d,y3d,z3d]=extrude3d(edgess,elements,tri,p(:,1),p(:,2),z0,z1,nz);
                p = [x3d,y3d,z3d];
            end
            
            %% Extrude upwards
%             dl=(z1-z0)/nze;
%             z = z0:dl:z1;
%             nnod = size(p,1);
%             znod = reshape( ones(nnod,1)*z, [], 1);
%             xynod = repmat(p,nz+1,1);
%             p = [xynod,znod];
%             mt = max(t(:));
%             
%             nzele = size(t,1);
%             prism = zeros(nzele*nz,6); % preallocate space for prism elements
%             l = (0:nz)*nzele;
%             c = 1;
%             for iz = 0:nz-1
%                 lo = l(c)+1;
%                 up = l(c+1);
%                 tt = [t(:,:)+mt*iz,t(:,:)+mt*(iz+1)];
%                 prism(lo:up,:) = tt;
%                 c=c+1;
%             end
% 
%             prism = unique(prism,'rows');
%             
%             %% Prism to tet
%             tet = [prism(:,[1,2,3,5]);
%                    prism(:,[1,5,6,4]); 
%                    prism(:,[1,5,6,3])]; 
            
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
                                        eli([4,1]);
                                        eli([1,3]);
                                        eli([2,4]);];

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
    
    methods(Hidden)
      function lh = addlistener(varargin)
         lh = addlistener@handle(varargin{:});
      end
      function notify(varargin)
         notify@handle(varargin{:});
      end
      function delete(varargin)
         delete@handle(varargin{:});
      end
      function Hmatch = findobj(varargin)
         Hmatch = findobj@handle(varargin{:});
      end
      function p = findprop(varargin)
         p = findprop@handle(varargin{:});
      end
      function TF = eq(varargin)
         TF = eq@handle(varargin{:});
      end
      function TF = ne(varargin)
         TF = ne@handle(varargin{:});
      end
      function TF = lt(varargin)
         TF = lt@handle(varargin{:});
      end
      function TF = le(varargin)
         TF = le@handle(varargin{:});
      end
      function TF = gt(varargin)
         TF = gt@handle(varargin{:});
      end
      function TF = ge(varargin)
         TF = ge@handle(varargin{:});
      end
   end
    
end


function [elements,edges,nodes]=tritoel(xnod,ynod,tri)
%columns of elements: longest edge, next (counterclock),next,next,
% oppnode, hangnode
%columns of edges: node1, node2, boundary (0 is interior edge),length

x2=xnod(tri(:,2))-xnod(tri(:,1));
y2=ynod(tri(:,2))-ynod(tri(:,1));
x3=xnod(tri(:,3))-xnod(tri(:,1));
y3=ynod(tri(:,3))-ynod(tri(:,1));
detj=x2.*y3-x3.*y2; %Area of elements
%%% Makes sure that the order is counterclockwise
ind=find(detj<0);
t=tri(ind,3);
tri(ind,3)=tri(ind,2);tri(ind,2)=t;
%%% 
nele=length(tri(:))/3; %number of elements

%% Find edges
edges=sort([tri(:,1),tri(:,2)],2);
edges=[edges;sort([tri(:,2),tri(:,3)],2)];
edges=[edges;sort([tri(:,1),tri(:,3)],2)];
edges=[edges,zeros(length(edges(:,1)),1),zeros(length(edges(:,1)),1)];

% Tr = TriRep(tri,xnod,ynod);
% edges2 = Tr.edges

%% Find doublets and remove them, add 1 or zero to the 3rd coulumn depending
%  on if that row did not contain double
%  Add length of edge to the fourth column

i1=1;
while i1 <= length(edges) %Size of edges changes...
	ind=find(edges(i1,1)==edges(:,1) & edges(i1,2)==edges(:,2));
	if(length(ind)==1)
		edges(i1,3)=1;
		n1=edges(i1,1);n2=edges(i1,2);
		edges(i1,4)=sqrt((xnod(n1)-xnod(n2))^2+(ynod(n1)-ynod(n2))^2);
	else
		n1=edges(ind(1),1);n2=edges(ind(1),2);
		edges(ind(1),4)=sqrt((xnod(n1)-xnod(n2))^2+(ynod(n1)-ynod(n2))^2);
		ind(1)=[];
		edges(ind,:)=[];
    end
    i1 = i1 + 1;

end
% edges
elements=zeros(nele,8);
for iel=1:nele
	ed1=sort([tri(iel,1),tri(iel,2)],2);
	ed1=find(ed1(1)==edges(:,1) & ed1(2)==edges(:,2));
	ed2=sort([tri(iel,3),tri(iel,2)],2);
	ed2=find(ed2(1)==edges(:,1) & ed2(2)==edges(:,2));
	ed3=sort([tri(iel,1),tri(iel,3)],2);
	ed3=find(ed3(1)==edges(:,1) & ed3(2)==edges(:,2));
	ed=[ed1,ed2,ed3,ed1,ed2];
	len=[edges(ed1,4),edges(ed2,4),edges(ed3,4)];
	ind=find(len==max(len));ind=ind(1);
	elements(iel,1:3)=[ed(ind),ed(ind+1),ed(ind+2)];
	nod=unique([edges(ed1,1:2),edges(ed2,1:2),edges(ed3,1:2)]);
	elements(iel,5)=setdiff(nod,edges(ed(ind),1:2));
end
nodes=[xnod,ynod];
end
function tri=trinode(edges,elements)
    %
    [nele,s]=size(elements);
    tri=zeros(nele,3);
    %el=elements(:,1:3);
    ed1=elements(:,1);ed2=elements(:,2);ed3=elements(:,3);
    e1=edges(ed1,1:2);e2=edges(ed2,1:2);e3=edges(ed3,1:2);
    % find number on e2 ~= number on e1
    ind=find(e1(:,1)~=e2(:,1) & e1(:,1)~=e2(:,2));
    tri1=e1(:,2);tri1(ind)=e1(ind,1);
    ind=find(e2(:,1)~=e3(:,1) & e2(:,1)~=e3(:,2));
    tri2=e2(:,2);tri2(ind)=e2(ind,1);
    ind=find(e3(:,1)~=e1(:,1) & e3(:,1)~=e1(:,2));
    tri3=e3(:,2);tri3(ind)=e3(ind,1);
    tri=[tri1,tri2,tri3];
end
function edges=direction(edges,xnod,ynod)
    nedges=size(edges,1);
    edges=edges(:,1:2);
    for ed=1:nedges
        v=[xnod(edges(ed,2))-xnod(edges(ed,1)),ynod(edges(ed,2))-ynod(edges(ed,1))];
        if(dot(v,[1,0]) < 0)
            edges(ed,1:2)=fliplr(edges(ed,1:2));
        end
    end
end
function [tri3d,x3d,y3d,z3d]=extrude3d(edges,elements,tri,xnod,ynod,zmin,zmax,nlayer)
    %
    %  4 ------ 6
    %   |\    /|
    %   | \  / |
    %   |  \/ 5|
    %   |   |  |
    %   |   |  |
    %   \   |  / 3
    %  1 \  | /
    %     \ |/
    %      \/
    %        2
    %  z=0 to z=1, built upwards
    %
    dl=(zmax-zmin)/nlayer;
    [nrtri,s]=size(tri);
    nno=length(xnod);
    tri3d=zeros(nrtri,4);
    x3d=zeros(2*nno,1);y3d=x3d;z3d=x3d;
    nele3=1;
    nod=0;
    for iel=1:nrtri
        edel=elements(iel,1:3);

        if(edges(edel(3),1)== edges(edel(1),1))
            nod=1;
        elseif(edges(edel(1),1) == edges(edel(2),1))
            nod=2;
        elseif(edges(edel(2),1)==edges(edel(3),1))
            nod=3;
        end
        switch nod
            case 1
                tri3d(nele3,1:4)=[tri(iel,1),tri(iel,2)+nno,tri(iel,3)+nno,tri(iel,1)+nno];
                if(edges(edel(2),1)== tri(iel,2))
                    tri3d(nele3+1,1:4)=[tri(iel,1:3),tri(iel,3)+nno];
                    tri3d(nele3+2,1:4)=[tri(iel,1),tri(iel,2),tri(iel,3)+nno,tri(iel,2)+nno]; % OK
                else
                    tri3d(nele3+1,1:4)=[tri(iel,1:3),tri(iel,2)+nno];
                    tri3d(nele3+2,1:4)=[tri(iel,1),tri(iel,2)+nno,tri(iel,3),tri(iel,3)+nno]; % OK
                end
                nele3=nele3+3;
            case 2
                tri3d(nele3,1:4)=[tri(iel,2),tri(iel,2)+nno,tri(iel,3)+nno,tri(iel,1)+nno];
                if(edges(edel(3),1)== tri(iel,3))
                    tri3d(nele3+1,1:4)=[tri(iel,1:3),tri(iel,1)+nno];
                    tri3d(nele3+2,1:4)=[tri(iel,1)+nno,tri(iel,2),tri(iel,3),tri(iel,3)+nno]; %OK
                else
                    tri3d(nele3+1,1:4)=[tri(iel,1:3),tri(iel,3)+nno];
                    tri3d(nele3+2,1:4)=[tri(iel,1)+nno,tri(iel,1),tri(iel,3)+nno,tri(iel,2)]; %OK
                end
                nele3=nele3+3;
            case 3
                tri3d(nele3,1:4)=[tri(iel,1)+nno,tri(iel,2)+nno,tri(iel,3),tri(iel,3)+nno];
                if(edges(edel(1),1)== tri(iel,1))
                    tri3d(nele3+1,1:4)=[tri(iel,1:3),tri(iel,2)+nno];
                    tri3d(nele3+2,1:4)=[tri(iel,1)+nno,tri(iel,2)+nno,tri(iel,1),tri(iel,3)]; %OK
                else
                    tri3d(nele3+1,1:4)=[tri(iel,1:3),tri(iel,1)+nno];
                    tri3d(nele3+2,1:4)=[tri(iel,1)+nno,tri(iel,2),tri(iel,3),tri(iel,2)+nno];
                end
                nele3=nele3+3;
        end
    end
    x3d(1:nno)=xnod;x3d(nno+1:end)=xnod;
    y3d(1:nno)=ynod;y3d(nno+1:end)=ynod;
    % z3d(1:nno)=0; %Original
    z3d(1:nno)=zmin;
    z3d(nno+1:end)=zmin+dl;
    z3dstart = unique(z3d);
    z3dstart = z3dstart(2);
    if(nlayer > 1)
            triex0=tri3d;

        for ilayer=2:nlayer
            triex=tri3d(:)+(ilayer-1)*nno;
            triex0=[triex0;reshape(triex,size(tri3d))];    
            x3d=[x3d;xnod];y3d=[y3d;ynod];
    %         z3d=[z3d;ones(size(xnod))*ilayer*dl]; % Original
            z3d=[z3d;ones(size(xnod))*(z3dstart+(ilayer-1)*dl)];
        end
        tri3d=triex0;
    end
end
