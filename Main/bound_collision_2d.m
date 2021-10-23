function collision_bound=bound_collision_2d(convp,domain_type,domain_size,domain_size2)
% 
collision_bound=0;
if strcmp(domain_type,'Circle')
    for i=1:size(convp,1)
        if (convp(i,1)^2+convp(i,2)^2)>domain_size^2
            collision_bound=1;
            break;
        end
    end
else
    for i=1:size(convp,1)
        if convp(i,1)<-domain_size/2 || convp(i,1)>domain_size/2 || convp(i,2)>domain_size2/2 || convp(i,2)<-domain_size2/2
            collision_bound=1;
            break;
        end
    end
end


end

