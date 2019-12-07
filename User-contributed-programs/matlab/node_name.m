%plot grid, to understand it
% M is the number of nodes,
% N is the number of triangles
% Jamie Pringle, University of New Hampshire

clf
clear all
load mesh

%plot nodes

col=0.9*[1 1 1];
jnk=patch('Vertices',mesh.nodexy/1e3,'Faces',mesh.trinodes,...
     'FaceColor',col,'EdgeColor',col);


for m=1:M
   text(mesh.nodexy(m,1)/1e3,mesh.nodexy(m,2)/1e3,sprintf('%4.0f',m),...
       'fontsize',6,'color','r');
end

% $$$ 
% $$$ for m=1:N
% $$$   text(mesh.uvnode(m,1)/1e3,mesh.uvnode(m,2)/1e3,sprintf('%4.0f',m),...
% $$$       'fontsize',6,'color','g');
% $$$ end

load depthgrd_fine

hold on

plot([mesh.nodexy(mesh.trinodes(:,1),1)' ...
      ;mesh.nodexy(mesh.trinodes(:,2),1)']/1e3,...
    [mesh.nodexy(mesh.trinodes(:,1),2)' ...
      ;mesh.nodexy(mesh.trinodes(:,2),2)']/1e3,'k')

plot([mesh.nodexy(mesh.trinodes(:,2),1)' ...
      ;mesh.nodexy(mesh.trinodes(:,3),1)']/1e3,...
    [mesh.nodexy(mesh.trinodes(:,2),2)' ...
      ;mesh.nodexy(mesh.trinodes(:,3),2)']/1e3,'k')

plot([mesh.nodexy(mesh.trinodes(:,3),1)' ...
      ;mesh.nodexy(mesh.trinodes(:,1),1)']/1e3,...
    [mesh.nodexy(mesh.trinodes(:,3),2)' ...
      ;mesh.nodexy(mesh.trinodes(:,1),2)']/1e3,'k')


%if true, plot depth, and pick out certain contours
if (1==2)
  [han,jnk]=contour(depthgrd.xvec/1e3,depthgrd.yvec/1e3,...
		    depthgrd.dpth,-[-270:20:0],'b-');
  [han,jnk]=contour(depthgrd.xvec/1e3,depthgrd.yvec/1e3,...
		    depthgrd.dpth,-[-140 -140],'r-');
 set(jnk,'linewidth',2)
 [han,jnk]=contour(depthgrd.xvec/1e3,depthgrd.yvec/1e3,...
		   depthgrd.dpth,-[-200 -200],'g-');
 set(jnk,'linewidth',2)
 for hand=jnk
   set(hand,'linewidth',2)
 end
end

%axis([1.0337    1.7947   -0.3238    0.3721]*1e3)

hold off
 
hold off
orient landscape

