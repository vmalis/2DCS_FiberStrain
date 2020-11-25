%=========================================================================
%   Step3 for 2d Strain Rate Analysis Toolbox (master script)
%
%       part of 2d Strain Rate Toolkit
%=========================================================================
%
%
% INput: (output of Step2_ROIs)
%
%           1)  results.mat
%           2)  study folder path
%
% OUTput:   1) plots in a folder 'results' created in a root folder of a
%              study. Default is pdf, eps option is commented.
%_____________________________________________________
% required subroutines:
%       1)rdir
%_____________________________________________________
% 
% v.2.0 tailored for ULLS data analysis
% (improved from v.1.0 used for Age-related Differences in Strain Rate Tensor paper (MRM 2014))
% 
% written by Vadim Malis
% 02/15 at UCSD RIL
%==========================================================================
%
% More info in HOWTO_2dSR.pdf
%

path=uigetdir();
cd(path)
mkdir('plots')

saving_path=[path '/plots'];

file_list=rdir('**/results.mat');

multiWaitbar('Creating plots', 0, 'Color', 'g');
N=size(file_list,1);


for i=1:N

    %load results file
    load(file_list(i).name)
    %number of frames
    n=size(NEV,2);               
    %geting all id paramters
    strtext=file_list(i).name;
    j=strfind(strtext,'_');
    k=strfind(strtext,'/');
    study_label=[strtext(1:min(j-1)) '-' strtext(j+1:min(k)-1)];
    subject_name=strtext(min(k)+1:min(k)+2);

    %Plots
    
    %---------------------------PEV----------------------------------------
    h_fig=figure('Visible','off');
    
    haxis1 = axes('Position',[0 0 1 1],'Visible','off');
    haxis2 = axes('Position',[0.07 0.2 0.9 .75]);
    
    intPEV1 = interp1(1:n,PEV(1,:)',1:0.1:n,'spline');
    intPEV2 = interp1(1:n,PEV(2,:)',1:0.1:n,'spline');
    intPEV3 = interp1(1:n,PEV(3,:)',1:0.1:n,'spline');
    
    h_p1=plot(PEV(1,:),'b.');
    hold on
    h_p2=plot(1:0.1:n,intPEV1,'b-');
    h_p3=plot(PEV(2,:),'r.');
    h_p4=plot(1:0.1:n,intPEV2,'r-');
    h_p5=plot(PEV(3,:),'g.');
    h_p6=plot(1:0.1:n,intPEV3,'g-');
    hold off
    legend([h_p2,h_p4,h_p6],'proximal', 'mid', 'distal');
    
    ylim([min(PEV(:)-100),max(PEV(:))+100]);
    xlim([0,n+1]);
    title('Positive eigen values')
    set(h_fig, 'Units', 'normalized', 'PaperPositionMode', 'auto');
           
    str(1) = {study_label};
    str(2) = {['Subject Initials: ' subject_name]};
    str(3) = {['Subject ID: ' PatientID]};
    str(3) = {['Slice Location: ' num2str(SliceLocation)]};
    str(4) = {['Series name: ' Series_name]};
    str(5) = {['Image Quality: ' quality]};
    set(gcf,'CurrentAxes',haxis1)
    text(.3,0.1,str,'FontSize',12)
    
    filename=[study_label '_' subject_name '_PEV_' num2str(int8(SliceLocation)) '.pdf'];
    fullname=[path '/plots/' filename];
%     print(h_fig, '-depsc2', '-r300', '-tiff', '-loose', fullname);
    print(h_fig, '-dpdf', '-r300', fullname);
    
    %---------------------------NEV----------------------------------------
    h_fig=figure('Visible','off');
    
    haxis1 = axes('Position',[0 0 1 1],'Visible','off');
    haxis2 = axes('Position',[0.07 0.2 0.9 .75]);
    
    intNEV1 = interp1(1:n,NEV(1,:),1:0.1:n,'spline');
    intNEV2 = interp1(1:n,NEV(2,:),1:0.1:n,'spline');
    intNEV3 = interp1(1:n,NEV(3,:),1:0.1:n,'spline');
    
    h_p1=plot(NEV(1,:),'b.');
    hold on
    h_p2=plot(1:0.1:n,intNEV1,'b-');
    h_p3=plot(NEV(2,:),'r.');
    h_p4=plot(1:0.1:n,intNEV2,'r-');
    h_p5=plot(NEV(3,:),'g.');
    h_p6=plot(1:0.1:n,intNEV3,'g-');
    hold off
    legend([h_p2,h_p4,h_p6],'proximal', 'mid', 'distal');
    
    ylim([min(NEV(:)-100),max(NEV(:))+100]);
    xlim([0,n+1]);
    title('Negative eigen values')
    set(h_fig, 'Units', 'normalized', 'PaperPositionMode', 'auto');
           
    str(1) = {study_label};
    str(2) = {['Subject Initials: ' subject_name]};
    str(3) = {['Subject ID: ' PatientID]};
    str(3) = {['Slice Location: ' num2str(SliceLocation)]};
    str(4) = {['Series name: ' Series_name]};
    str(5) = {['Image Quality: ' quality]};
    set(gcf,'CurrentAxes',haxis1)
    text(.3,0.1,str,'FontSize',12)
    
    filename=[study_label '_' subject_name '_NEV_' num2str(int8(SliceLocation)) '.pdf'];
    fullname=[path '/plots/' filename];
%     print(h_fig, '-depsc2', '-r300', '-tiff', '-loose', fullname);
    print(h_fig, '-dpdf', '-r300', fullname);

    %---------------------------SEV----------------------------------------
    h_fig=figure('Visible','off');
    
    haxis1 = axes('Position',[0 0 1 1],'Visible','off');
    haxis2 = axes('Position',[0.07 0.2 0.9 .75]);
    
    intSEV1 = interp1(1:n,SEV(1,:),1:0.1:n,'spline');
    intSEV2 = interp1(1:n,SEV(2,:),1:0.1:n,'spline');
    intSEV3 = interp1(1:n,SEV(3,:),1:0.1:n,'spline');
    
    h_p1=plot(SEV(1,:),'b.');
    hold on
    h_p2=plot(1:0.1:n,intSEV1,'b-');
    h_p3=plot(SEV(2,:),'r.');
    h_p4=plot(1:0.1:n,intSEV2,'r-');
    h_p5=plot(SEV(3,:),'g.');
    h_p6=plot(1:0.1:n,intSEV3,'g-');
    hold off
    legend([h_p2,h_p4,h_p6],'proximal', 'mid', 'distal');
    
    ylim([min(SEV(:)-100),max(SEV(:))+100]);
    xlim([0,n+1]);
    title('Sum eigen values')
    set(h_fig, 'Units', 'normalized', 'PaperPositionMode', 'auto');
           
    str(1) = {study_label};
    str(2) = {['Subject Initials: ' subject_name]};
    str(3) = {['Subject ID: ' PatientID]};
    str(3) = {['Slice Location: ' num2str(SliceLocation)]};
    str(4) = {['Series name: ' Series_name]};
    str(5) = {['Image Quality: ' quality]};
    set(gcf,'CurrentAxes',haxis1)
    text(.3,0.1,str,'FontSize',12)
    
    filename=[study_label '_' subject_name '_SEV_' num2str(int8(SliceLocation)) '.pdf'];
    fullname=[path '/plots/' filename];
%     print(h_fig, '-depsc2', '-r300', '-tiff', '-loose', fullname);
    print(h_fig, '-dpdf', '-r300',fullname);
    %---------------------------SR ANGL------------------------------------
    h_fig=figure('Visible','off');
    
    haxis1 = axes('Position',[0 0 1 1],'Visible','off');
    haxis2 = axes('Position',[0.07 0.2 0.9 .75]);
    
    intANGL1 = interp1(1:n,ANGL(1,:),1:0.1:n,'spline');
    intANGL2 = interp1(1:n,ANGL(2,:),1:0.1:n,'spline');
    intANGL3 = interp1(1:n,ANGL(3,:),1:0.1:n,'spline');
    
    h_p1=plot(ANGL(1,:),'b.');
    hold on
    h_p2=plot(1:0.1:n,intANGL1,'b-');
    h_p3=plot(ANGL(2,:),'r.');
    h_p4=plot(1:0.1:n,intANGL2,'r-');
    h_p5=plot(ANGL(3,:),'g.');
    h_p6=plot(1:0.1:n,intANGL3,'g-');
    hold off
    legend([h_p2,h_p4,h_p6],'proximal', 'mid', 'distal');
    
    ylim([-100,100]);
    xlim([0,n+1]);
    title('SR Angle')
    set(h_fig, 'Units', 'normalized', 'PaperPositionMode', 'auto');
           
    str(1) = {study_label};
    str(2) = {['Subject Initials: ' subject_name]};
    str(3) = {['Subject ID: ' PatientID]};
    str(3) = {['Slice Location: ' num2str(SliceLocation)]};
    str(4) = {['Series name: ' Series_name]};
    str(5) = {['Image Quality: ' quality]};
    set(gcf,'CurrentAxes',haxis1)
    text(.3,0.1,str,'FontSize',12)
    
    filename=[study_label '_' subject_name '_SR-ANG_' num2str(int8(SliceLocation)) '.pdf'];
    fullname=[path '/plots/' filename];
%     print(h_fig, '-depsc2', '-r300', '-tiff', '-loose', fullname);
    print(h_fig, '-dpdf', '-r300',fullname);

    %---------------------------FIBER ANGL---------------------------------
    h_fig=figure('Visible','off');
    
    haxis1 = axes('Position',[0 0 1 1],'Visible','off');
    haxis2 = axes('Position',[0.07 0.2 0.9 .75]);
    
    intFiber_ANGL1 = interp1(1:n,Fiber_ANGL(1,:),1:0.1:n,'spline');
    intFiber_ANGL2 = interp1(1:n,Fiber_ANGL(2,:),1:0.1:n,'spline');
    intFiber_ANGL3 = interp1(1:n,Fiber_ANGL(3,:),1:0.1:n,'spline');
    
    h_p1=plot(Fiber_ANGL(1,:),'b.');
    hold on
    h_p2=plot(1:0.1:n,intFiber_ANGL1,'b-');
    h_p3=plot(Fiber_ANGL(2,:),'r.');
    h_p4=plot(1:0.1:n,intFiber_ANGL2,'r-');
    h_p5=plot(Fiber_ANGL(3,:),'g.');
    h_p6=plot(1:0.1:n,intFiber_ANGL3,'g-');
    hold off
    legend([h_p2,h_p4,h_p6],'proximal', 'mid', 'distal');
    
    ylim([min(Fiber_ANGL(:))-5,max(Fiber_ANGL(:))+5]);
    xlim([0,n+1]);
    title('Fiber Angle')
    set(h_fig, 'Units', 'normalized', 'PaperPositionMode', 'auto');
           
    str(1) = {study_label};
    str(2) = {['Subject Initials: ' subject_name]};
    str(3) = {['Subject ID: ' PatientID]};
    str(3) = {['Slice Location: ' num2str(SliceLocation)]};
    str(4) = {['Series name: ' Series_name]};
    str(5) = {['Image Quality: ' quality]};
    set(gcf,'CurrentAxes',haxis1)
    text(.3,0.1,str,'FontSize',12)
    
    filename=[study_label '_' subject_name '_Fiber-ANG_' num2str(int8(SliceLocation)) '.pdf'];
    fullname=[path '/plots/' filename];
%     print(h_fig, '-depsc2', '-r300', '-tiff', '-loose', fullname);
    print(h_fig, '-dpdf', '-r300', fullname);
    
multiWaitbar('Creating plots', i/N);

end

multiWaitbar('Creating plots', 'Close');  
