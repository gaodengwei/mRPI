function h = sector(x0,y0,u1,v1,u2,v2,r)
plot(x0,y0,'k.')
hold on
theta1 = cart2pol(u1,v1);
theta2 = cart2pol(u2,v2);
h = [];

for i=1:size(x0,1)
    for j=1:size(x0,2)
        t = linspace(theta1(i,j),theta2(i,j),100);
        if theta1(i,j)-theta2(i,j)>pi 
            t = linspace(theta1(i,j),theta2(i,j)+2*pi,100);
%         t=t;
        end
        x = x0(i,j)+r*cos(t);
        y = y0(i,j)+r*sin(t);
        h0 = fill([x0(i,j),x,x0(i,j)],[y0(i,j),y,y0(i,j)],'b','EdgeColor','none');
        set(h0,'FaceColo', [0.5 0.5 0.5],'FaceAlpha',0.2)
        h = [h,h0];
        
    end
end
end