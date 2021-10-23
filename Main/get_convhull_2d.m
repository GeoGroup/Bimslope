function convp=get_convhull_2d(p)
%
k = convhull(p(:,1),p(:,2));
convp=[p(k,1),p(k,2)];
convp(length(convp),:)=[];
end

