function [edges,nedges]=insert_edge(pfrom,pto,edges,nedges)
% 
% modify edges to connect node pfrom to pto
%
% edges a (M,6) matrix, containing the 6 nodes each node is connected to
% nedeges a (M) vector, with the number of edges each node is already connected
%                       to.
%
% lousy, slow algorythem, but easy to code, and assumes little about inputs
% Jamie Pringle, University of New Hampshire

if (nedges(pfrom)~=6)
  for ne = 1:nedges(pfrom)
    if (edges(pfrom,ne)==pto)
      return %edge already in, return
    end
  end
  edges(pfrom,nedges(pfrom)+1)=pto;
  nedges(pfrom)=nedges(pfrom)+1;
else
  %need not do anything, found all edges
  return
end

