function pos_all=update_pos_2d(poly,pos_all)
% delete the points in the polyhedron 
% index=in_polyhedron(K,P,pos_all); % this is slow
index=inpolygon(pos_all(:,1),pos_all(:,2),poly(:,1),poly(:,2));
% pos_allt=pos_all(index,:);
% plot(pos_allt(:,1),pos_allt(:,2),'r.');hold on
pos_all=pos_all(~index,:);

end

