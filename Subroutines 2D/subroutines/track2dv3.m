function [xs,ys,vr,vx,vy,vz] = track2dv3(x,y,v_x,v_y,v_z,dt,res)
%   -----------------------------------------------------------------------
%   usage: [xs,ys,vr,vx,vy,vz] = track2dv3(x,y,v_x,v_y,v_z,dt,res)
%   -----------------------------------------------------------------------
%   remarks: track2dv3 tracks tissue displacement in 2d using pc-mri data
%   -----------------------------------------------------------------------
%   input arguments:
%   ----------------
%   x,y: seed point coordinates (num_seeds x 1)
%   v_x: calibrated velocity in horizontal direciton
%   v_y: calibrated velocity in vertical direction
%   v_z: calibrated velocity in through-plane direction
%   dt: time difference b/w consecutive temporal phases
%   res: image pixel resolution [pixels/cm]
%   -----------------------------------------------------------------------
%   returns:
%   --------
%   xs,ys: coordinates of trajectories (num_seeds, x num_phases)
%   vr: resultant velocity of trajectories (num_seeds, x num_phases)
%   vx: horizontal component velocity (num_seeds, x num_phases)
%   vy: vertical component velocity (num_seeds, x num_phases)
%   vz: through-plane component velocity (num_seeds, x num_phases)
%   -----------------------------------------------------------------------
%   David Shin
%   09/07/2008 - v3
%       - Removed 3x3 velocity averaging. 
%       - Use averaging of v_x, v_y, and v_z prior to running this
%         function instead. 
%       - Added 23th phase in xs and ys. 
%   -----------------------------------------------------------------------
%   Vadim Malis
%   02/07/2015  - v4
%       -Catch for points going out of image boundaries
%       -Waitbar is added
%

%% pre-allocate variables
numphases = size(v_y,3);
xs = zeros(size(x,1),numphases);
ys = zeros(size(y,1),numphases);
vr = zeros(size(y,1),numphases);
vx = zeros(size(y,1),numphases);
vy = zeros(size(y,1),numphases);
vz = zeros(size(y,1),numphases);

%% averaging v_x, v_y, v_z
h = fspecial('average', [3,3]);
for i = 1:numphases;
    v_x(:,:,i) = imfilter(v_x(:,:,i),h);
    v_y(:,:,i) = imfilter(v_y(:,:,i),h);
    v_z(:,:,i) = imfilter(v_z(:,:,i),h);
end

%% row column dims
row=size(v_x,1);
column=size(v_x,2);

multiWaitbar('tracking...', 0, 'Color', 'g');

%% track engine
for i = 1:numphases
    

    if (i == 1)
        xs(:,i) = x(:);
        ys(:,i) = y(:);
    else
        for j = 1:length(x)
            slice = i - 1;
            
            if isnan(xs(j,i-1))
                xs(j,i-1)=xs(j,i-2);
            end
            
            if isnan(ys(j,i-1))
                ys(j,i-1)=ys(j,i-2);
            end
            
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

            
            z_x(1,1) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), slice);
            z_x(2,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), slice);
            z_x(1,2) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), slice);
            z_x(2,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), slice);
            
            z_y(1,1) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), slice);
            z_y(2,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), slice);
            z_y(1,2) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), slice);
            z_y(2,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), slice);

            z_z(1,1) = v_z(uint16(y_temp_floor), uint16(x_temp_floor), slice);
            z_z(2,1) = v_z(uint16(y_temp_ceil), uint16(x_temp_floor), slice);
            z_z(1,2) = v_z(uint16(y_temp_floor), uint16(x_temp_ceil), slice);
            z_z(2,2) = v_z(uint16(y_temp_ceil), uint16(x_temp_ceil), slice);
            %-------------------------------------------------------
            v_x_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_x, px, py,'linear');
            v_y_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_y, px, py,'linear');
            v_z_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_z, px, py,'linear');
            %-------------------------------------------------------
            delta_x = v_x_temp * dt(i-1) * res;
            delta_y = v_y_temp * dt(i-1) * res;

            xs(j,i) = xs(j,i-1) - delta_x;
            ys(j,i) = ys(j,i-1) + delta_y;
            

            
            vx(j,i-1) = v_x_temp;
            vy(j,i-1) = v_y_temp;
            vz(j,i-1) = v_z_temp;
            vr(j,i-1) = sqrt(v_x_temp^2 + v_y_temp^2);

            if (i == numphases) %this segment is needed to find velocity at the last frame
                
                if isnan(xs(j,numphases))
                xs(j,numphases)=xs(j,numphases-1);
                end
            
                if isnan(ys(j,numphases))
                ys(j,numphases)=ys(j,numphases-1);
                end
                
                px = xs(j,numphases);
                py = ys(j,numphases);

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
                
                if x_temp_floor == x_temp_ceil
                    x_temp_ceil = x_temp_floor + 1;
                end
                if y_temp_floor == y_temp_ceil
                    y_temp_ceil = y_temp_floor + 1;
                end
                
                z_x(1,1) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), numphases);
                z_x(2,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), numphases);
                z_x(1,2) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), numphases);
                z_x(2,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), numphases);

                z_y(1,1) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), numphases);
                z_y(2,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), numphases);
                z_y(1,2) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), numphases);
                z_y(2,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), numphases);

                z_z(1,1) = v_z(uint16(y_temp_floor), uint16(x_temp_floor), slice);
                z_z(2,1) = v_z(uint16(y_temp_ceil), uint16(x_temp_floor), slice);
                z_z(1,2) = v_z(uint16(y_temp_floor), uint16(x_temp_ceil), slice);
                z_z(2,2) = v_z(uint16(y_temp_ceil), uint16(x_temp_ceil), slice);
                %-------------------------------------------------------
                v_x_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_x, px, py,'linear');
                v_y_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_y, px, py,'linear');
                v_z_temp = interp2([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil], z_z, px, py,'linear');
                %-------------------------------------------------------
                vx(j,numphases) = v_x_temp;
                vy(j,numphases) = v_y_temp;
                vz(j,numphases) = v_z_temp;
                vr(j,numphases) = sqrt(v_x_temp^2 + v_y_temp^2);
            end % last phase velocity calculation
        end % j loop (# of seeds)
    end % if
    
    
    multiWaitbar('tracking...', i/numphases);
    
end % i loop (# of phases)

multiWaitbar('tracking...', 'Close');  

