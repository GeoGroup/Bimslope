function [p] = particle_generate_polygon_2d(rmin,rmax,dr,eta,n1,n2,ext,rot)
% Generate a single 2D particle
% r0=10;
% dr=1;
% eta=0.7;
% n=15;
r0=rmin+(rmax-rmin)*rand;
dr=r0*dr;
n=randi([n1,n2],1,1);
% if rmin<0.75
%     n=6;
% end
r=zeros(n,1);
dtheta=zeros(n,1);
for i=1:n
    sj=rand(1,1);
    r(i)=(r0+(2*sj-1)*dr);
    dtheta(i)=2*pi/n+(2*sj-1)*eta*2*pi/n;
end
dtheta=dtheta*2*pi/sum(dtheta);
dtheta=cumsum(dtheta);

% r=r.*(1+(ext-1)*abs(cos(dtheta)));
x=r.*cos(dtheta); %+rot/pi
y=r.*sin(dtheta); %+rot/pi
x=x*ext;
% r=sqrt(x.^2+y.^2);
% x=r.*cos(dtheta); %+rot/pi
% y=r.*sin(dtheta); %+rot/pi
if rot>=0
    xn=x*cos(-rot/180*pi)-y*sin(-rot/180*pi);
    yn=x*sin(-rot/180*pi)+y*cos(-rot/180*pi);
else
    rot=rand(1,1)*pi;
    xn=x*cos(-rot)-y*sin(-rot);
    yn=x*sin(-rot)+y*cos(-rot);
end
p=[xn,yn];
end