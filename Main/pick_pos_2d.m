function pos=pick_pos_2d(pos_all)
% pick position from all points
% n=unidrnd(size(pos_all,1));
% pos=pos_all(n,:);
pos=pos_all(unidrnd(size(pos_all,1)),:);
end

