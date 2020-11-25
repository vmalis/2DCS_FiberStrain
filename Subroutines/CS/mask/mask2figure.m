function []=mask2figure(K_pattern,file_name_img)

% Generates undersampling pattern
%   
% Input:   
%      kspace
%
% Output:  
%      mask
%
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis

% save undersampling pattern to image
figure
imshow(repelem(K_pattern',1,4),'InitialMagnification',800)
axis on
ax = gca;
xticks(1:4:size(K_pattern,1)*4);
    
labels = string(1:size(K_pattern,1));
n=1;

    for i=1:4:size(labels,1)
        labels(i)=string(n);
        n=n+1;
    end
    
ax.XAxis.TickLabels = labels; % set
set(gca,'TickLabelInterpreter','latex');
xlabel('frame number','Interpreter','latex');
ylabel('phase encoding','Interpreter','latex');
set(gcf,'color','w');
box off
export_fig(file_name_img)
close
end