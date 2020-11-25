




%crop image along x
t=find(squeeze(abs(sum(test_m(:,:,1),1))));
if isempty(t)
x1=1;
x2=size(m,2);
else
x1=t(1);
x2=t(end);
end
mag = test_m(:,x1:x2,:);



figure("Position",[10,10,2000,400],'color','k')
ax1 = axes;
M=reshape(mag,[size(mag,1),size(mag,2)*size(mag,3)]);
imagesc(M)
colormap(ax1,'gray');
daspect([1 1 1])
axis image
h=gca;
h.Visible="off";

ax2 = axes;
imagesc(ax2,flip(M,2),'alphadata',0.3);
cmap=jet(100);
cmap(1:10,:) = 0;
colormap(ax2,cmap);
caxis(ax2,[0 10]);
h=gca;
h.Visible="off";
axis image

daspect([1 1 1])
linkprop([ax1 ax2],'Position');
h=colorbar;
ylabel(h, '$E_{octahedral} \quad \mathrm{[mm]}$','fontsize',16,'interpreter','latex','color','white')
h.Color = 'w';
h.TickLabelInterpreter='latex';
h.Limits = [0,10];

export_fig figure1.png -m16