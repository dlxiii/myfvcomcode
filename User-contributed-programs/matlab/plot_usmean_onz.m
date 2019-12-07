% plot salinity and currents on a z-surface
% Jamie Pringle, University of New Hampshire

clf
clear all

%location of data files
filedir='data/'
filehist='chn_history.nc';

zplot=-5; %depth to plot
nt=1;      %time to plot (in records of output in file) 
dnt=1;     %number of output records to average over 
ax=1e5*[0.8302 5.4016 0 1.9548];     %domain to plot, km. An empty matrix plots
	   %the entire domain
	   
cax=[29 30.04]; %color range of plot

dts=8.64e4*6.0; %the number of seconds so that the length of the
                %current arrows = dts*U 

min_plot=0.01;  %do not plot arrows where current is less than
                %min_plot
		
scale_arrow_len = 0.1; %length of scale current arrow in m/s		

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

%interp salinity to depth surface
 tic;
 sdepth=nan*ones(size(squeeze(schunk(1,:))));
 for m=1:M
   %if (rem(m,1000)==0)
   %  m
   %end
   sdepth(m)=interp1_fast(mesh.nodez(:,m),schunk(:,m),zplot);
 end
 toc
 
%interp u and v to depth, but only on points to be ploted. 
load minspace1000  %this is a file of points to be ploted from trim_uv.m
gp=goodpts;
u=zeros(length(gp),1);
v=zeros(length(gp),1);
tic
for g=1:length(gp)
  %if (rem(g,1000)==0)
  %  g
  %end
  u(g)=interp1_fast(mesh.zuv(:,gp(g)),uchunk(:,gp(g)),zplot);
  v(g)=interp1_fast(mesh.zuv(:,gp(g)),vchunk(:,gp(g)),zplot);
end
toc

%draw coastlines
 col=0.7*[1 1 1];
 jnk=patch('Vertices',mesh.nodexy,'Faces',mesh.trinodes,...
     'FaceColor',col,'EdgeColor',col);

%plot salinity
 patch('Vertices',mesh.nodexy,'Faces',mesh.trinodes,'Cdata',sdepth,...
     'edgecolor','interp','facecolor','interp')

 
 suptitle(sprintf(...
     'U and S at z=%2.0f for file %s at time %4.1f hours',...
     zplot,filehist,time))

 if ~isempty(ax)
   axis(ax)
 end
 
 caxis(cax)
 axis equal
 
%draw arrows

 %load what points to draw arrows on
 %load allspace
 gp=goodpts;
 
 dx=dts*u';
 dy=dts*v';
 offset=[dx;dy]';

 %don't plot arrows where their are no currents, or where they are less
 %than min_plot
 if ~isempty(ax)
   plotscale=1;
   bp=find((~isnan(u)).*(min_plot<sqrt(u.^2+v.^2)).*...
	   (mesh.uvnode(gp,1)>ax(1)*plotscale).*...
	   (mesh.uvnode(gp,1)<ax(2)*plotscale).*...
	   (mesh.uvnode(gp,2)>ax(3)*plotscale).*...
	   (mesh.uvnode(gp,2)<ax(4)*plotscale));
 else
   bp=find((~isnan(u)).*(min_plot<sqrt(u.^2+v.^2)));
 end
 
 if ~isempty(ax)
   axis(ax)
 end
 
 arrow(mesh.uvnode(gp(bp),:),(mesh.uvnode(gp(bp),:)+offset(bp,:)),...
     'length',4,'tipangle',20,'width',0.5);
 
 %plot a scale arrow.
 ax_now=axis;
 xa=ax_now(1)+0.1*(ax_now(2)-ax_now(1));
 ya=ax_now(3)+0.9*(ax_now(4)-ax_now(3));
 ya_txt=ax_now(3)+0.93*(ax_now(4)-ax_now(3));
 
 arrow([xa ya],[xa+dts*scale_arrow_len ya],...
       'length',4,'tipangle',20,'width',0.5);
 text(xa,ya_txt,sprintf('%4.2f cm/s',scale_arrow_len*100)) 
 
 
 
 title(sprintf(['Arrow is distance traveled in %2.1f days,'...
		' min arrow speed=%4.3fm/s'],dts/8.64e4,min_plot))

 if ~isempty(ax)
   axis(ax)
 end
 
 %draw depth contours if the test below is true
 if (1==2)
   load depthgrd
   hold on
   contour(depthgrd.xvec,depthgrd.yvec,...
	   depthgrd.dpth,[0:5:100],'k-')
   hold off
   if ~isempty(ax)
     axis(ax)
   end
 end
 
%You need this for some buggy versions of matlab, with some buggy
%OpenGL drivers.  
 set(gcf,'renderer','painters')

colorbar


%clean up
close(nchist)

