%load model mesh in 
%
% Jamie Pringle, University of New Hampshire
%
% Saves out several .mat files with geometric information
% 
% mesh.mat contains
%   M the number of nodes, on which all scalers are stored (T,S,eta, H etc)
%   N the number of triangles, roughly twice M
%     also the number of U,V points, each on the centroid of the triangle
%   KBM1 the number of sigma levels on which U,V,T,S etc are stored.
%   
%   and the structure "mesh", whose feilds are
%   mesh.surround(N,3) the U,V (or centroid) points clockwise arround 
%                      a U,V point, called NBE in model and paper
%   mesh.trinodes(N,3) the t/s nodes surrounding a triangle, in clockwise order
%                      called NV in model and N(j) in paper
%   mesh.edges(M,6)  the up to 6 other nodes that each node is connected
%                    to. NaN if the edge does not exist.
%   mesh.nedges(M,1) the number of edges each node is connected to.
%   mesh.nodexy(M,2) the x and y position of the nodes, x is 1 coloumn
%   mesh.depth(M) the depth
%   mesh.uvnode(N,2) the x and y positions of the centroids of the triangles
%                    at which U and V are defined. (sorry for the poor name)
%   mesh.nodez(KBM1,M) the depth of each scaler, assuming free surface=0
%   mesh.zuv(KBM1,N)   the depth of each U,V point, assuming free surface=0
%   mesh.sigvec(KBM1) the sigma coordinate grid
%
% interpcoef.mat contains
% the interpolation coefficients to interpolate u,v and the scalers
% from where they are stored to any surrounding point. see the routines 
% window_sect_inf_make.f and window_sect_data_make.f for information, or
% secflux.m for how they are used.
%
% depthgrd.mat contains depth gridded onto a course rectangular grid
%    depthgrd.xvec is x dimension of grid
%    depthgrd.yvec is y dimension of grid
%    depthgrid.dpth is the depth
%
% depthgrd_fine.mat contains depth gridded onto a finer rectangular grid
%
%

clear all

%filename of datafile and grid file
 filegrid='data/chn_setup.nc';

%get M, the number of nodes, and N, the number of triangles
 ncgrid=netcdf(filegrid,'nowrite');  
 M=length(ncgrid{'VX'}); %number of nodes
 N=length(ncgrid{'XC'}); %number of triangles
 KBM1=size(ncgrid{'ZALL'},1);
 
%get mesh data
 mesh.surround=ncgrid{'NBE'}(:)';
 mesh.trinodes=ncgrid{'NV'}(1:3,:)';

%read in node positions
 mesh.nodexy=zeros(M,2);
 mesh.nodexy(:,1)=ncgrid{'VX'}';
 mesh.nodexy(:,2)=ncgrid{'VY'}';

%get bathy data, and adjust as in model
 mesh.depth=ncgrid{'H'}(:);

%form u,v locations (the centroid of triangle)
 uvnodes=zeros(N,2);
 uvdepth=zeros(N,1);
 for n=1:N
   uvnode(n,1)=(mesh.nodexy(mesh.trinodes(n,1),1)+...
       mesh.nodexy(mesh.trinodes(n,2),1)+...
       mesh.nodexy(mesh.trinodes(n,3),1))/3;
   uvnode(n,2)=(mesh.nodexy(mesh.trinodes(n,1),2)+...
       mesh.nodexy(mesh.trinodes(n,2),2)+...
       mesh.nodexy(mesh.trinodes(n,3),2))/3;
   uvdepth(n)=(mesh.depth(mesh.trinodes(n,1))+...
       mesh.depth(mesh.trinodes(n,2))+...
       mesh.depth(mesh.trinodes(n,3)))/3;
 end
 mesh.uvnode=uvnode;
 mesh.uvdepth=uvdepth;

%get data on how each node is connected to three other nodes
 edges=zeros(M,6)*nan;
 nedges=zeros(M,1); %number of edges already found
 for n=1:N
   n;
   [edges,nedges]=insert_edge(mesh.trinodes(n,1),mesh.trinodes(n,2),...
       edges,nedges);
   [edges,nedges]=insert_edge(mesh.trinodes(n,2),mesh.trinodes(n,3),...
       edges,nedges);
   [edges,nedges]=insert_edge(mesh.trinodes(n,3),mesh.trinodes(n,1),...
       edges,nedges);
   [edges,nedges]=insert_edge(mesh.trinodes(n,2),mesh.trinodes(n,1),...
       edges,nedges);
   [edges,nedges]=insert_edge(mesh.trinodes(n,3),mesh.trinodes(n,2),...
       edges,nedges);
   [edges,nedges]=insert_edge(mesh.trinodes(n,1),mesh.trinodes(n,3),...
       edges,nedges);
 end
 mesh.edges=edges;
 mesh.nedges=nedges;

%create a full xyz list of node locations
%
%for node, their are KBM1 depths, so we must include all these depths in ...
%the list of locations
 zall=repmat((1:KBM1)',1,M);
 zuv=repmat((1:KBM1)',1,N);

%make z in depth coordinates, not the discrete equivalent of sigman coords.
 dsig=1/(KBM1+1);
 sigvec=(0:(KBM1-1))'*dsig+0.5*dsig;
 for m=1:M
   zall(:,m)=-sigvec*mesh.depth(m);
 end
 for n=1:N
   zuv(:,n)=-sigvec*mesh.uvdepth(n);
 end
 mesh.nodez=zall;
 mesh.zuv=zuv;
 mesh.sigvec=sigvec;
 
 save mesh mesh M N KBM1

    
%the following produces depth on a 2D regular grid
%for ease of visualization.

%define the regular grid
 ndx=100;
 xvec=min(mesh.nodexy(:,1)):(max(mesh.nodexy(:,1))- ...
			     min(mesh.nodexy(:,1)))/ndx:max(mesh.nodexy(:,1)); 
 yvec=min(mesh.nodexy(:,2)):(max(mesh.nodexy(:,2))- ...
			     min(mesh.nodexy(:,2)))/ndx:max(mesh.nodexy(:,2)); ;
 
%make the 2D depth array, and save
 dpth=griddata(mesh.nodexy(:,1),mesh.nodexy(:,2),mesh.depth,xvec,yvec');
 depthgrd.xvec=xvec;
 depthgrd.yvec=yvec;
 depthgrd.dpth=dpth;
 save depthgrd depthgrd

%define the regular grid
 ndx=300;
 xvec=min(mesh.nodexy(:,1)):(max(mesh.nodexy(:,1))- ...
			     min(mesh.nodexy(:,1)))/ndx:max(mesh.nodexy(:,1)); 
 yvec=min(mesh.nodexy(:,2)):(max(mesh.nodexy(:,2))- ...
			     min(mesh.nodexy(:,2)))/ndx:max(mesh.nodexy(:,2)); ;
 
%make the 2D depth array, and save
 dpth=griddata(mesh.nodexy(:,1),mesh.nodexy(:,2),mesh.depth,xvec,yvec');
 depthgrd.xvec=xvec;
 depthgrd.yvec=yvec;
 depthgrd.dpth=dpth;
 save depthgrd_fine depthgrd


 
clf
%plot nodes
plot(mesh.nodexy(:,1),mesh.nodexy(:,2),'b.',...
     mesh.uvnode(:,1),mesh.uvnode(:,2),'r.');

hold on
for m=1:M
  for ne=1:nedges(m)
    plot([mesh.nodexy(m,1) mesh.nodexy(edges(m,ne),1)],...
	[mesh.nodexy(m,2) mesh.nodexy(edges(m,ne),2)],'k-')
  end
end
hold off
    

