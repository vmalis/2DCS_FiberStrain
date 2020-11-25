function [] = coilmap2figure(CM,raw,n)

% Saves coilmap to image
%   
% Input:   
%      CM  - coilmap
%      raw - pfilestruct
%
% Output:  
%      -
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis

CM_images=zeros(size(CM,1),size(CM,1),size(CM,3));

file_list=rdir('*.7');
p=GERecon('Pfile.Load',file_list(n).name);
GERecon('Pfile.SetActive', p)

    for coil=1:size(CM,3)
        % correct for gradient nonlinearity
        temp=gradwarp(CM(:,:,coil),raw.corners(1));
        CM_images(:,:,coil) = abs(temp);
    end
   
    %prepare for crop
    t=find(squeeze(abs(CM_images(end,:,1))));
    if isempty(t)
        x1=1;
        x2=size(CM_images,2);
    else
        x1=t(1);
        x2=t(end);
    end
    Images=zeros(size(CM,1),x2-x1+1,3,size(CM,3));
    
    for coil=1:size(CM,3)
        %lable and crop
        lable{1}=['coil #' num2str(coil)];
        Images(:,:,:,coil) = insertText(abs(CM_images(:,x1:x2,coil)),[0,200],lable,'FontSize',18,'BoxColor',...
    'green','BoxOpacity',0.4,'TextColor','white');
        
    end

    montage(Images,'Size',[2,ceil(size(Images,4)/2)],'BackgroundColor','white','BorderSize',[2,2])
    hF = gcf;
    hF.Position(3:4) = [size(Images,2)*size(Images,4)/2 size(Images,1)*2+60];
    
    title('Normalized coil sensitivity','Interpreter','latex','FontSize',32,'Color','white');
    set(gcf,'color','black');
    colormap('gray')
    h=colorbar;
    h.Color='w';
    h.FontSize=20;
    set(h,'TickLabelInterpreter','latex');
    % clearvars p temp lable CM_images
    filename=[raw.series(1:3),'coilmap.png'];
    export_fig(filename,'-m8')
    close all
    
    
end