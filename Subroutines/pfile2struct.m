function [pFileSt]=pfile2struct


%PathName = uigetdir('~/Desktop','choose p-files directory');
%cd(PathName)
file_list=rdir('*.7')

for i=1:size(file_list,1)
    i
    pinfo=GERecon('Pfile.Load',file_list(i).name);
    
    GERecon('Pfile.SetActive',pinfo)
    header=GERecon('Pfile.Header');
        
    pFileSt(i).filename=pinfo.filename;
    pFileSt(i).series=header.SeriesData.se_desc;
    pFileSt(i).header=header;
    pFileSt(i).slices=pinfo.slices;
    pFileSt(i).xres=pinfo.xRes;
    pFileSt(i).yres=pinfo.yRes;
    pFileSt(i).frames=pinfo.phases;
    pFileSt(i).echoes=pinfo.echoes;
    pFileSt(i).channels=pinfo.channels;
    
    
    
    for k=1:pFileSt(i).slices
        
        pFileSt(i).corners(k) = GERecon('Pfile.Corners', k);
        
        for l=1:pFileSt(i).echoes
            for m=1:pFileSt(i).channels
                for n=1:pFileSt(i).frames
        
                pFileSt(i).kspace(:,:,k,l,m,n) = GERecon('Pfile.KSpace',k,l,m,n);
    
                end
            end
        end
    end
    
    
end

