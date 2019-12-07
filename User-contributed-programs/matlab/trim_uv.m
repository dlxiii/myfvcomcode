% this code saves the indices of u,v points at least minspace apart from each
% other, to keep the graphics from getting cluttered. 
%
%
% saves the following files
% allspace.mat is all u/v points in a file
% minspace*.mat is a collection of u/v points at least minspace
% apart
%
% Jamie Pringle, University of New Hampshire


minspace=100000; % meters

load mesh

goodpts=1:N;
save allspace goodpts

for n=1:N
  if ~isnan(goodpts(n))
    dist=sqrt((mesh.uvnode(:,1)-mesh.uvnode(n,1)).^2+...
	(mesh.uvnode(:,2)-mesh.uvnode(n,2)).^2);
    indx=find(dist<minspace);
    goodpts(indx)=nan;
    goodpts(n)=n;
    sum(~isnan(goodpts));
    n
  end
end
total_left=sum(~isnan(goodpts))

goodpts=goodpts(find(~isnan(goodpts)));

save(sprintf('minspace%0.1d',minspace),'goodpts');


plot(mesh.uvnode(:,1)/1e3,mesh.uvnode(:,2)/1e3,'r.',...
    mesh.uvnode(goodpts,1)/1e3,mesh.uvnode(goodpts,2)/1e3,'b.');  
