function [x,temp_average,v]=tempPCA(data)

FltAr=array_fltn3d(data);
temp_average = mean(FltAr,2);
FltAr_subt=FltAr-temp_average;
FltAr_subt_T=ctranspose(FltAr_subt);
C=double(FltAr_subt_T*FltAr_subt);


[v,~]=eigs(C,size(C,1));
u=FltAr_subt*v;
x=reshape(u,size(data));

end