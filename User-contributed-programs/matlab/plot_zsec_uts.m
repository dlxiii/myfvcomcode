% plot vertical cross-section of U/T/S
% Jamie Pringle, University of New Hampshire
clear all
clf

%location of data files
filedir='data/'
filehist='chn_history.nc';

zplot=-5; %depth to make (x,y) plot
nt=1;      %time to plot (in records of output in file) 
dnt=1;     %number of output records to average over 
ax=1e5*[0.8302 5.4016 0 1.9548];     %domain to plot in (x,y),
	       %km. An empty matrix plots 
	       %the entire domain
	   
cax_s=[29 30.04]; %color range of salt
cax_t=[15 25]; %color range of temp
cax_u=[-20:1:20]; %contour intervals of velocity in CM/S!

dsec=1e3; %this code intepolates the data onto a line in x,y
          %space.  dsec is the spacing between points on the line,
          %and should be much less than the node spacing for most
          %accurate results. 


%load grid data
load mesh

%attach netcdf files to variables
nchist=netcdf([filedir filehist],'nowrite');
t=nchist{'t1'};
s=nchist{'s1'};
u=nchist{'u'};
v=nchist{'v'};
el=nchist{'el'};
thour=nchist{'thour'};


%get data
tchunk=squeeze(t(nt:(nt+dnt-1),:,:));
schunk=squeeze(s(nt:(nt+dnt-1),:,:));
uchunk=squeeze(u(nt:(nt+dnt-1),:,:));
vchunk=squeeze(v(nt:(nt+dnt-1),:,:));
time=thour(nt:(nt+dnt-1));

%take time mean from nt:(nt+dnt-1)
if (dnt~=1)
  time=squeeze(mean(time));
  tchunk=squeeze(mean(tchunk,1));
  schunk=squeeze(mean(schunk,1));
  uchunk=squeeze(mean(uchunk,1));
  vchunk=squeeze(mean(vchunk,1));
end

%interp salinity to depth in order to plot horizontal map
 tic;
 sdepth=nan*ones(size(squeeze(schunk(1,:))));
 for m=1:M
   sdepth(m)=interp1_fast(mesh.nodez(:,m),schunk(:,m),zplot);
 end
 toc

%plot map with salinity

subplot(2,1,1)

%plot vertices
 col=0.7*[1 1 1];
 jnk=patch('Vertices',mesh.nodexy/1,'Faces',mesh.trinodes,...
	   'FaceColor',col,'EdgeColor',col);
 
 hold on
 patch('Vertices',mesh.nodexy/1,'Faces',mesh.trinodes,'Cdata',sdepth,...
     'edgecolor','interp','facecolor','interp')
 hold off
 caxis(cax_s);

 
 title(sprintf('S at z=%2.0f for run %s at time %4.1f hours',...
     zplot,filehist,time))
 
 
%add depth to plot
 if (1==2)
   load depthgrd
   hold on
   contour(depthgrd.xvec/1,depthgrd.yvec/1,...
	   depthgrd.dpth,[0:10:100],'k-')
   hold off
 end

 axis equal
 if ~isempty(ax)
   axis(ax)
 end
 colorbar

%save axis
 map=gca;

%get two points to define end points of section
 [xp,yp]=ginput(2)

 % xp=[0.8720    1.0112]*1e6
 % yp=[3.9117   -9.8075]*1e4
 
%define section to calculate along
 dg=dsec/sqrt((xp(2)-xp(1)).^2+(yp(2)-yp(1)).^2);
 gvec=(0:dg:1)';
 xvec=(xp(2)-xp(1))*gvec+xp(1);
 yvec=(yp(2)-yp(1))*gvec+yp(1);
 distvec=sqrt((xvec-xp(1)).^2+(yvec-yp(1)).^2);

 %compute angle of line, and angle perpendicular to line
 theta=atan2(yp(2)-yp(1),xp(2)-xp(1));
 cross_theta=theta-pi/2;
 
 hold on
 plot(xvec/1,yvec/1,'b*',xvec/1,yvec/1,'y-')
 hold off
 
%interpolate T,S, U and V and Z onto section on k level.
 tic;
 echo on
 tsec=griddata_vect(mesh.nodexy(:,1),mesh.nodexy(:,2),...
     tchunk,xvec,yvec);
 ssec=griddata_vect(mesh.nodexy(:,1),mesh.nodexy(:,2),...
     schunk,xvec,yvec);
 zsec=griddata_vect(mesh.nodexy(:,1),mesh.nodexy(:,2),...
     mesh.nodez,xvec,yvec);
 usec=griddata_vect(mesh.uvnode(:,1),mesh.uvnode(:,2),...
     uchunk,xvec,yvec);
 vsec=griddata_vect(mesh.uvnode(:,1),mesh.uvnode(:,2),...
     vchunk,xvec,yvec);
 echo off
 toc 

 %make matrix of horizontal distances for use in makeing plots
 distsec=[];
 for k=1:KBM1
   distsec(k,:)=distvec';
 end

 %make cross section velocity.
 %positive velocity is out of section
 cross_vel=sin(cross_theta)*vsec+cos(cross_theta)*usec;
 
 subplot(2,2,3)
 pcolor(distsec/1,zsec,tsec);
 hold on; 
 [clab,hlab]=contour(distsec/1,zsec,cross_vel*100,[cax_u],'w-');
 clabel(clab,hlab);
 negpick(clab,hlab,'r','w','k');
 hold off
 shading flat
 xlabel('temperature')
 ylabel('depth')
 caxis(cax_t)
 colorbar
 shading interp
 
 title('contours are velocity, positive is out of page');

 
 drawnow
 
 subplot(2,2,4)
 pcolor(distsec/1,zsec,ssec);
 shading flat
 hold on; 
 [clab,hlab]=contour(distsec/1,zsec,cross_vel*100,[cax_u],'w-');
 clabel(clab,hlab);
 negpick(clab,hlab,'r','w','k');
 hold off
 xlabel('salinity')
 ylabel('depth')
 caxis(cax_s)
 colorbar
 
 shading interp

 title('contours are velocity, positive is out of page');
 
 
 suptitle(sprintf(...
     'tidal mean values for run %s at time %4.1f hours',...
     filehist,time))

 orient landscape
 
 %You need this for some buggy versions of matlab, with some buggy
 %OpenGL drivers.  
 set(gcf,'renderer','painters')

 