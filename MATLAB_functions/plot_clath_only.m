%C is a matrix with coordinates of clathrin. 
%fnum is the Figure number.
%Assumes COM is the first row.
function[sigma, n1]=plot_clath_only(fnum, C)
size(C)

figure(fnum);
hold on;
cmap=colormap;
cleg=C(2,:)-C(1,:);
LegLen=sqrt(cleg*cleg')
full=zeros(4,3);
for i=2:1:4
    half=C(i,:)-C(1,:);
    full(i,:)=half*2.0+C(1,:);
    
    p=[C(1,:); C(i,:)];
    
    size(p)
    plot3(p(:,1),p(:,2),p(:,3),'d-','Color',cmap(i*2+0,:),'LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');%CLA binding sites
end
for i=5:1:7
    p=[C(1,:); C(i,:)];
  plot3(p(:,1),p(:,2),p(:,3),'o-','Color',cmap(i*3+30,:),'LineWidth',4);%AP2 binding sites
end

