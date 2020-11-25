function []=vepc_error_plot(dVx,dVy,dVz,name_suffix)

%montage of velocities

lbl_space=zeros(size(dVx,1),size(dVx,2));

I_montage=cat(3,cat(3,lbl_space,dVx),cat(3,lbl_space,dVy),cat(3,lbl_space,dVz));
cmap=jet(1000);
cmap(1,:) = 0;
border_w=20;
montage(I_montage,'DisplayRange',[1 50],'Size',[3,size(dVx,3)+1],'BorderSize',[border_w,0],'ThumbnailSize',[size(dVx,1),size(dVx,2)]);
colormap(cmap)
set(gcf,'color','black');
h=colorbar;
h.Color='w';
h.FontSize=12;
ylabel(h, '[RMSE %]','Interpreter','latex','FontSize',20);
set(h,'TickLabelInterpreter','latex','FontSize',20); 

text(round(size(dVx,2)/3),round(border_w+size(dVx,1)/2),'$v_{x}$','FontSize',32,'Interpreter','latex','Color','white')
text(round(size(dVx,2)/3),round(3*border_w+size(dVx,1)*1.5),'$v_{y}$','FontSize',32,'Interpreter','latex','Color','white')
text(round(size(dVx,2)/3),round(5*border_w+size(dVx,1)*2.5),'$v_{z}$','FontSize',32,'Interpreter','latex','Color','white')

f=msgbox('Saving to image...', '','help');
filename=[name_suffix,'.png'];
export_fig(filename,'-m8')

close(f)
close all