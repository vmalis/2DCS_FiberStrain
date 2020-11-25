function [peak, sd] = peak_strains(varargin)

    x             =   varargin{1};
    SD            =   varargin{2};

    if nargin==2

        [~,idx]    =   max(abs(varargin{1}));
            if size(size(x),2)==3
                for i=1:size(x,2)
                    for j=1:size(x,3)   
                        peak(i,j)=x(idx(1,i,j),i,j);
                        sd(i,j)=SD(idx(1,i,j),i,j);
                    end
                end
                
            elseif size(size(x),2)==2
                for i=1:size(x,2)   
                        peak(i)=x(idx(1,i),i);
                        sd(i)=SD(idx(1,i),i);
                end
            else
                peak=x(idx);
                sd=SD(idx);
            end
        
        
    else
    
        frame = varargin{3};
        
        if size(size(x),2)==3
                for i=1:size(x,2)
                    for j=1:size(x,3)   
                        peak(i,j)=x(frame,i,j);
                        sd(i,j)=SD(frame,i,j);
                    end
                end
                
            elseif size(size(x),2)==2
                for i=1:size(x,2)   
                        peak(i)=x(frame,i);
                        sd(i)=SD(frame,i);
                end
            else
                peak=x(frame);
                sd=SD(frame);
         end
        
    end


end