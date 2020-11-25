SR_ratio_iso_post=[SR_data(1:51).SR_ratio];
SR_ratio_iso_post=reshape(SR_ratio_iso_post,[51,22,3]);

for i=1:size(SR_ratio_iso_post,2)
    for j=1:size(SR_ratio_iso_post,3)
        A=squeeze(SR_ratio_iso_post(:,i,j));
        SR_ratio_iso_post_filt(i,j)=stdfilter(A);
    end
end
stdev=squeeze(nanstd(SR_ratio_iso_post,1));
plots_4paper(SR_ratio_iso_post_filt,stdev,'SR_f_f/SR_f_c','SR_f_f/SR_f_c','SR_ratio_iso_post_abs')






SR_ratio_iso_post=[SR_data(65:112).SR_ratio];
SR_ratio_iso_post=reshape(SR_ratio_iso_post,[48,22,3]);

for i=1:size(SR_ratio_iso_post,2)
    for j=1:size(SR_ratio_iso_post,3)
        A=squeeze(SR_ratio_iso_post(:,i,j));
        SR_ratio_iso_post_filt(i,j)=stdfilter(A);
    end
end
stdev=squeeze(nanstd(SR_ratio_iso_post,1));
plots_4paper(SR_ratio_iso_post_filt,stdev,'SR_f_f/SR_f_c','SR_f_f/SR_f_c','SR_ratio_iso_pre_abs')





SR_ratio_iso_post=[SR_data(1:52).SR_ratio];
SR_ratio_iso_post=reshape(SR_ratio_iso_post,[52,22,3]);

for i=1:size(SR_ratio_iso_post,2)
    for j=1:size(SR_ratio_iso_post,3)
        A=squeeze(SR_ratio_iso_post(:,i,j));
        SR_ratio_iso_post_filt(i,j)=stdfilter(A);
    end
end
stdev=squeeze(nanstd(SR_ratio_iso_post,1));
plots_4paper(SR_ratio_iso_post_filt,stdev,'SR_f_f/SR_f_c','SR_f_f/SR_f_c','SR_ratio_ecc_post_abs')






SR_ratio_iso_post=[SR_data(68:116).SR_ratio];
SR_ratio_iso_post=reshape(SR_ratio_iso_post,[49,22,3]);

for i=1:size(SR_ratio_iso_post,2)
    for j=1:size(SR_ratio_iso_post,3)
        A=squeeze(SR_ratio_iso_post(:,i,j));
        SR_ratio_iso_post_filt(i,j)=stdfilter(A);
    end
end
stdev=squeeze(nanstd(SR_ratio_iso_post,1));
plots_4paper(SR_ratio_iso_post_filt,stdev,'SR_f_f/SR_f_c','SR_f_f/SR_f_c','SR_ratio_ecc_pre_abs')