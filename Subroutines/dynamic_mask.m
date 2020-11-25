function [mask_dynamic] = dynamic_mask(xs,ys,zs,mask)

mask(isnan(mask))=0;

xs=round(xs.*repmat(mask,[1,1,3,size(xs,4)]));
ys=round(ys.*repmat(mask,[1,1,3,size(ys,4)]));
zs=round(zs.*repmat(mask,[1,1,3,size(zs,4)]));

xs(xs>size(xs,2))=size(xs,2);
ys(ys>size(ys,1))=size(ys,1);
zs(zs>size(zs,3))=size(zs,3);

xs(xs<1)=1;
ys(ys<1)=1;
zs(zs<1)=1;
% 
% size(xs)
% size(ys)
% size(zs)

mask_dynamic=zeros(size(xs));

for t=1:size(xs,4)
    for y=1:size(xs,1)
        for x=1:size(xs,2)
            for z=1:size(xs,3)
                
                if xs(y,x,z,t)~=0 && ys(y,x,z,t)~=0 && zs(y,x,z,t)~=0
                    
                    mask_dynamic(ys(y,x,z,t),xs(y,x,z,t),zs(y,x,z,t),t) = 1;
                
                end
                
                
            end
        end
    end 
end




end