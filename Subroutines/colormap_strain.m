function []=colormap_strain(m,d,suffix,limits,units)


%crop image along x
t=find(squeeze(abs(sum(m(:,:,1),1))));
if isempty(t)
    x1=1;
    x2=size(m,2);
else
    x1=t(1);
    x2=t(end);
end
mag = m(:,x1:x2,:);
data= d(:,x1:x2,:);

%flatten for montage
M=reshape(mag,[size(mag,1),size(mag,2)*size(mag,3)]);
D=reshape(data,[size(data,1),size(data,2)*size(data,3)]);

figure("Position",[10,10,2000,400],'color','k')

ax1 = axes;
imagesc(M)
colormap(ax1,'gray');
daspect([1 1 1])
h=gca;
h.Visible="off";

ax2 = axes;
imagesc(ax2,D,'alphadata',0.2);

cmap=jet(301);

    if limits(1)==0
        cmap(1,:) = 0;
    elseif limits(1)<0 && limits(2)>0
        cmap(151,:) = 0;
    elseif limits(1)<0 && limits(2)==0
        cmap(end,:) =0;
    end    


colormap(ax2,cmap);
caxis(ax2,[limits(1) limits(2)]);
h=gca;
h.Visible="off";
daspect([1 1 1])
linkprop([ax1 ax2],'Position');

h=colorbar;

ylabel(h, units,'fontsize',16,'interpreter','latex','color','white')
h.Color = 'w';
h.TickLabelInterpreter='latex';
h.Limits = [limits(1) limits(2)];

filename=[suffix,'.png'];
export_fig(filename, '-m16')

close all

end