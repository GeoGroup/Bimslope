function edge = edge_from_vertice(point,np,label)
%Generate edge from node
edge=zeros(size(point,1),3)+label;
for i=1:size(point,1)
    if i==size(point,1)
        edge(i,1:2)=[np+i,np+1];
    else
        edge(i,1:2)=[np+i,np+i+1];
    end
end
end

