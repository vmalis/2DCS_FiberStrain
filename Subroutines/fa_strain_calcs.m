function [ED,LD,SRD] = fa_strain_calcs(sr_struct,dti_vec)



ED=zeros(size(sr_struct.E));
LD=ED;
SRD=ED;


for y=1:size(sr_struct.E,1)

    for x=1:size(sr_struct.E,2)
        
        for z=1
            
            for t=1:size(sr_struct.E,4)
                
               ED(y,x,z,t,:,:)=dti_rotation(dti_vec(y,x,z,:,:),sr_struct.E(y,x,z,t,:,:));
               LD(y,x,z,t,:,:)=dti_rotation(dti_vec(y,x,z,:,:),sr_struct.L(y,x,z,t,:,:));
               SRD(y,x,z,t,:,:)=dti_rotation(dti_vec(y,x,z,:,:),sr_struct.SR(y,x,z,t,:,:));
                
                
            end
        end
    end
end


end