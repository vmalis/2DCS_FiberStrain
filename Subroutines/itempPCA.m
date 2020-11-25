function [X]=itempPCA(x,temp_average,v)

u=reshape(x,size(x,1)*size(x,2),size(x,3));

X=u*v'+temp_average;
X=reshape(X,size(x));

end