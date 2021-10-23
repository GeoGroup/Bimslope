function [ b,l,w] = get_bound_2d( p )
% Get the boundary of polygon
% (y-y2)/(y1-y2) = (x-x2)/(x1-x2)
% 1/(x1-x2)*x+[-1/(y1-y2)]*y+y2/(y1-y2)-x2/(x1-x2)
b=zeros(5,2);
n=size(p,1);
D=pdist(p);
D=squareform(D);
T=max(D(:));
[row,col]=find(D==T);
row=sort(row);
x1=p(row(1),1);x2=p(row(2),1);
y1=p(row(1),2);y2=p(row(2),2);
l=sqrt((x2-x1)^2+(y2-y1)^2);
A=1/(x1-x2);B=-1/(y1-y2);C=y2/(y1-y2)-x2/(x1-x2);
up=row(1)+1:row(2)-1;
all=setdiff(1:n,row);
down=setdiff(all,up);
dist=zeros(length(up),1);
maxv=0;maxi=1;
for i=1:length(up)
    dist(i)=(A*p(up(i),1)+B*p(up(i),2)+C)/sqrt(A^2+B^2);
    if abs(dist(i))>maxv
        maxv=abs(dist(i));
        maxi=i;
    end
end
w=maxv;
x1n=x1+dist(maxi)*A/sqrt(A^2+B^2);y1n=y1+dist(maxi)*B/sqrt(A^2+B^2);
x2n=x2+dist(maxi)*A/sqrt(A^2+B^2);y2n=y2+dist(maxi)*B/sqrt(A^2+B^2);
% hold on;plot([x1n,x2n],[y1n,y2n],'-*');
b(1,1)=x2n;b(1,2)=y2n;
b(2,1)=x1n;b(2,2)=y1n;

b(5,1)=b(1,1);b(5,2)=b(1,2);
maxv=0;maxi=1;
for i=1:length(down)
    dist(i)=(A*p(down(i),1)+B*p(down(i),2)+C)/sqrt(A^2+B^2);
    if abs(dist(i))>maxv
        maxv=abs(dist(i));
        maxi=i;
    end
end
w=w+maxv;
x1n=x1+dist(maxi)*A/sqrt(A^2+B^2);y1n=y1+dist(maxi)*B/sqrt(A^2+B^2);
x2n=x2+dist(maxi)*A/sqrt(A^2+B^2);y2n=y2+dist(maxi)*B/sqrt(A^2+B^2);
% hold on;plot([x1n,x2n],[y1n,y2n],'-*');
b(4,1)=x2n;b(4,2)=y2n;
b(3,1)=x1n;b(3,2)=y1n;

end

