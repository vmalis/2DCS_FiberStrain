function [xs,ys] = track2dv4(x,y,v_x,v_y,dt,res,start_frame)

%   -----------------------------------------------------------------------
%   usage: [xs,ys] = track2dv4(x,y,v_x,v_y,dt,res,start_frame)
%   -----------------------------------------------------------------------
%   remarks: track2dv4 tracks tissue displacement in 2d using pc-mri data
%   -----------------------------------------------------------------------
%   input arguments:
%   ----------------
%   x,y: seed point coordinates (num_seeds x 1)
%   x:   col
%   y:   row
%   v_x: calibrated velocity in horizontal direciton
%   v_y: calibrated velocity in vertical direction
%   v_z: calibrated velocity in through-plane direction
%   dt: time difference b/w consecutive temporal phases
%   res: image pixel resolution [pixels/cm]
%   start_frame: initial frame to start tracking
%   -----------------------------------------------------------------------
%   returns:
%   --------
%   xs,ys: coordinates of trajectories (num_seeds, x num_phases)
%   -----------------------------------------------------------------------
%   David Shin
%   09/07/2008 - v3
%       - Removed 3x3 velocity averaging. 
%       - Use averaging of v_x, v_y, and v_z prior to running this
%         function instead. 
%       - Added 23th phase in xs and ys. 
%   -----------------------------------------------------------------------
%   Vadim Malis
%   02/22/2015  - v4
%       -Catch for points going out of image boundaries
%       -Waitbar is added
%       -Start frame (script allows tracking from arbitary frame)
%
%% pre-allocate variables

numphases = size(v_y,3);
xs = zeros(size(x,1),numphases+1);
ys = zeros(size(y,1),numphases+1);


%% averaging v_x, v_y, v_z
h = fspecial('average', [3,3]);
for i = 1:numphases;
    v_x(:,:,i) = imfilter(v_x(:,:,i),h);
    v_y(:,:,i) = imfilter(v_y(:,:,i),h);
end

multiWaitbar('tracking...', 0, 'Color', 'g');

%% track engine
T=0;

xs(:,start_frame+1) = x(:);
ys(:,start_frame+1) = y(:);

for i=start_frame:-1:1

T=T+1;    

        for j = 1:length(x)
            frame = i+1;
            
            px = xs(j,i+1);
            py = ys(j,i+1);
            
            % subpixel interpolation
            x_temp_floor = floor(px);
            x_temp_ceil = ceil(px);
            y_temp_floor = floor(py);
            y_temp_ceil = ceil(py);
            
            
            if x_temp_floor<=1
               x_temp_floor=1;
               x_temp_ceil=2;
            end
            
            if y_temp_floor<=1
               y_temp_floor=1;
               y_temp_ceil=2;
            end
            
            
            
            % catch for floor and ceil to be differnt
            if (x_temp_floor == x_temp_ceil)
                x_temp_ceil = x_temp_floor + 1;
            end
            if (y_temp_floor == y_temp_ceil)
                y_temp_ceil = y_temp_floor + 1;
            end

            
            z_x(1,1) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), frame);
            z_x(2,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), frame);
            z_x(1,2) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), frame);
            z_x(2,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), frame);
            
            z_y(1,1) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), frame);
            z_y(2,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), frame);
            z_y(1,2) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), frame);
            z_y(2,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), frame);

            %-------------------------------------------------------
            v_x_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_x, px, py,'linear');
            v_y_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_y, px, py,'linear');
            %-------------------------------------------------------
            delta_x = v_x_temp * dt(i) * res;
            delta_y = v_y_temp * dt(i) * res;

            xs(j,i) = xs(j,i+1) + delta_x;
            ys(j,i) = ys(j,i+1) - delta_y;

        end % j loop (# of seeds)
    
    multiWaitbar('tracking...', T/numphases);
    
end % i loop (# of phases)

xs(:,1)=[];
ys(:,1)=[];    
    
if start_frame<numphases    
  
    
for i = start_frame+1:numphases
T=T+1;    
        for j = 1:length(x)
            frame = i-1;
            
            px = xs(j,i-1);
            py = ys(j,i-1);
            
            % subpixel interpolation
            x_temp_floor = floor(px);
            x_temp_ceil = ceil(px);
            y_temp_floor = floor(py);
            y_temp_ceil = ceil(py);
            
            
            if x_temp_floor<=1
               x_temp_floor=1;
               x_temp_ceil=2;
            end
            
            if y_temp_floor<=1
               y_temp_floor=1;
               y_temp_ceil=2;
            end
            

            % catch for floor and ceil to be differnt
            if (x_temp_floor == x_temp_ceil)
                x_temp_ceil = x_temp_floor + 1;
            end
            if (y_temp_floor == y_temp_ceil)
                y_temp_ceil = y_temp_floor + 1;
            end

            z_x(1,1) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), frame);
            z_x(2,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), frame);
            z_x(1,2) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), frame);
            z_x(2,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), frame);
            
            z_y(1,1) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), frame);
            z_y(2,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), frame);
            z_y(1,2) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), frame);
            z_y(2,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), frame);

            %-------------------------------------------------------
            v_x_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_x, px, py,'linear');
            v_y_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_y, px, py,'linear');
            %-------------------------------------------------------
            delta_x = v_x_temp * dt(i-1) * res;
            delta_y = v_y_temp * dt(i-1) * res;

            xs(j,i) = xs(j,i-1) - delta_x;
            ys(j,i) = ys(j,i-1) + delta_y;
            

        end % j loop (# of seeds)
        multiWaitbar('tracking...', T/numphases);
    
end % i loop (# of phases)
    
    
end 


multiWaitbar('tracking...', 'Close');  