function stats=CSMVC_univariate(Y,data,fields)

    for i=1:size(fields,2)
        stats(i).variable=fields{i};
        temp=fitlm([data.(fields{i})],Y);
        stats(i).r=sqrt(temp.Rsquared.Ordinary)*sign(temp.Coefficients.Estimate(2));
        stats(i).absr=sqrt(temp.Rsquared.Ordinary);
        stats(i).p=temp.Coefficients.pValue(2);
        
    end
    
stats=sortStruct(stats,'absr',-1);
    
end