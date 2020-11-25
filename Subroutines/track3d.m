function [xs,ys,zs,vr,vx,vy,vz] = track3d(x,y,z,v_x,v_y,v_z,dt,pxs)
%   -----------------------------------------------------------------------
%   usage: [xs,ys,vr,vx,vy,vz] = track3d(x,y,z,v_x,v_y,v_z,dt,pxs)
%   -----------------------------------------------------------------------
%   input arguments:
%   ----------------
%   x,y,z: seed point coordinates
%   v_x: calibrated velocity in horizontal direciton
%   v_y: calibrated velocity in vertical direction
%   v_z: calibrated velocity in through-plane direction
%   dt: time difference b/w consecutive temporal phases
%   pxs: pixel spacing (in code is converted to resolution = pixels/cm)
%   -----------------------------------------------------------------------
%   returns:
%   --------
%   xs,ys,zs: coordinates of trajectories
%   vr: velocity vector magnitude
%   vx: horizontal component velocity
%   vy: vertical component velocity
%   vz: through-plane component velocity
%   -----------------------------------------------------------------------
%   Extended version of tracking algorythm written by David Shin
%   _____________________________________________________
%   written by Vadim Malis
%   10/14 at UCSD RIL
%   -----------------------------------------------------------------------

%% pre-allocate variables
numphases = size(v_y,4);
xs = zeros(size(x,1),numphases);
ys = zeros(size(x,1),numphases);
zs = zeros(size(x,1),numphases);
vr = zeros(size(x,1),numphases);
vx = zeros(size(x,1),numphases);
vy = zeros(size(x,1),numphases);
vz = zeros(size(x,1),numphases);

multiWaitbar('Tracking...', 0, 'Color', 'g');
waitn=numphases*length(x);
step=0;


size_x=size(v_x,1);
size_y=size(v_x,2);
size_z=size(v_x,3);


%% track engine
for i = 1:numphases

    
    if (i == 1)
        xs(:,i) = x(:);
        ys(:,i) = y(:);
        zs(:,i) = z(:);
    else
        for j = 1:length(x) %through all points inside ROI
            
            step=step+1;
            multiWaitbar('Tracking...', step/waitn);
            
            
            slice = i - 1;
            
            %initial coordinates, coming from previous frame
            px = xs(j,i-1);
            py = ys(j,i-1);
            pz = zs(j,i-1);
            
            % subpixel interpolation
            x_temp_floor = floor(px);
            x_temp_ceil = ceil(px);
            y_temp_floor = floor(py);
            y_temp_ceil = ceil(py);
            z_temp_floor = floor(pz);
            z_temp_ceil = ceil(pz);

            % catch for floor and ceil to be differnt
            if (x_temp_floor == x_temp_ceil)
                x_temp_ceil = x_temp_floor + 1;
            end
            
            if (y_temp_floor == y_temp_ceil)
                y_temp_ceil = y_temp_floor + 1;
            end
            
            if (z_temp_floor == z_temp_ceil)
                z_temp_ceil = z_temp_floor + 1;
            end
            
            %out of bounds check
            if x_temp_floor<=1
               x_temp_floor=1;
               x_temp_ceil=2;
            end
            
            if x_temp_ceil>=size_x
               x_temp_ceil=size_x;
               x_temp_floor=size_x-1;
            end
            
            if y_temp_floor<=1
               y_temp_floor=1;
               y_temp_ceil=2;
            end
            
            if y_temp_ceil>=size_y
               y_temp_ceil=size_y;
               y_temp_floor=size_y-1;
            end
            
            
            if z_temp_floor<=1
               z_temp_floor=1;
               z_temp_ceil=2;
            end
            
            if z_temp_ceil>=size_z
               z_temp_ceil=size_z;
               z_temp_floor=size_z-1;
            end

            
            z_x(1,1,1) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_floor), slice);
            z_x(1,1,2) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_ceil), slice);
            z_x(1,2,1) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_floor), slice);
            z_x(1,2,2) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_ceil), slice);
            z_x(2,1,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_floor), slice);
            z_x(2,1,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_ceil), slice);
            z_x(2,2,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_floor), slice);
            z_x(2,2,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_ceil), slice);
            
            z_y(1,1,1) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_floor), slice);
            z_y(1,1,2) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_ceil), slice);
            z_y(1,2,1) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_floor), slice);
            z_y(1,2,2) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_ceil), slice);
            z_y(2,1,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_floor), slice);
            z_y(2,1,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_ceil), slice);
            z_y(2,2,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_floor), slice);
            z_y(2,2,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_ceil), slice);
            
            z_z(1,1,1) = v_z(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_floor), slice);
            z_z(1,1,2) = v_z(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_ceil), slice);
            z_z(1,2,1) = v_z(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_floor), slice);
            z_z(1,2,2) = v_z(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_ceil), slice);
            z_z(2,1,1) = v_z(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_floor), slice);
            z_z(2,1,2) = v_z(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_ceil), slice);
            z_z(2,2,1) = v_z(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_floor), slice);
            z_z(2,2,2) = v_z(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_ceil), slice);
            
            %-------------------------------------------------------
            v_x_temp = interp3([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil],[z_temp_floor z_temp_ceil], z_x, px, py, pz,'linear');
            v_y_temp = interp3([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil],[z_temp_floor z_temp_ceil], z_y, px, py, pz,'linear');
            v_z_temp = interp3([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil],[z_temp_floor z_temp_ceil], z_z, px, py, pz,'linear');
            
                            
                if isnan(v_x_temp)
                    v_x_temp=0;
                    
                end
                if isnan (v_y_temp)
                    v_y_temp=0;
                end
                if isnan (v_z_temp)
                    v_z_temp=0;
                end
            
            
            %-------------------------------------------------------
            delta_x = v_x_temp * dt(i-1) * 1/(pxs(1)/10);
            delta_y = v_y_temp * dt(i-1) * 1/(pxs(2)/10);
            delta_z = v_z_temp * dt(i-1) * 1/(pxs(3)/10);
            
            
            xs(j,i) = xs(j,i-1) + delta_x;
            ys(j,i) = ys(j,i-1) + delta_y;
            zs(j,i) = zs(j,i-1) + delta_z;
            
            vx(j,i-1) = v_x_temp;
            vy(j,i-1) = v_y_temp;
            vz(j,i-1) = v_z_temp;
            vr(j,i-1) = sqrt(v_x_temp^2 + v_y_temp^2+v_z_temp^2);

            if (i == numphases) %this segment is needed to find velocity at the last frame
                px = xs(j,numphases);
                py = ys(j,numphases);
                pz = zs(j,numphases);
                
                % subpixel interpolation
                x_temp_floor = floor(px);
                x_temp_ceil = ceil(px);
                y_temp_floor = floor(py);
                y_temp_ceil = ceil(py);
                z_temp_floor = floor(pz);
                z_temp_ceil = ceil(pz);
                
                
            % catch for floor and ceil to be differnt
            if (x_temp_floor == x_temp_ceil)
                x_temp_ceil = x_temp_floor + 1;
            end
            
            if (y_temp_floor == y_temp_ceil)
                y_temp_ceil = y_temp_floor + 1;
            end
            
            if (z_temp_floor == z_temp_ceil)
                z_temp_ceil = z_temp_floor + 1;
            end
            
            %out of bounds check
            if x_temp_floor<=1
               x_temp_floor=1;
               x_temp_ceil=2;
            end
            
            if x_temp_ceil>=size_x
               x_temp_ceil=size_x;
               x_temp_floor=size_x-1;
            end
            
            if y_temp_floor<=1
               y_temp_floor=1;
               y_temp_ceil=2;
            end
            
            if y_temp_ceil>=size_y
               y_temp_ceil=size_y;
               y_temp_floor=size_y-1;
            end
            
            
            if z_temp_floor<=1
               z_temp_floor=1;
               z_temp_ceil=2;
            end
            
            if z_temp_ceil>=size_z
               z_temp_ceil=size_z;
               z_temp_floor=size_z-1;
            end
                
                z_x(1,1,1) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_floor), numphases);
                z_x(1,1,2) = v_x(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_ceil), numphases);
                z_x(1,2,1) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_floor), numphases);
                z_x(1,2,2) = v_x(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_ceil), numphases);
                z_x(2,1,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_floor), numphases);
                z_x(2,1,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_ceil), numphases);
                z_x(2,2,1) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_floor), numphases);
                z_x(2,2,2) = v_x(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_ceil), numphases);
            
                z_y(1,1,1) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_floor), numphases);
                z_y(1,1,2) = v_y(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_ceil), numphases);
                z_y(1,2,1) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_floor), numphases);
                z_y(1,2,2) = v_y(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_ceil), numphases);
                z_y(2,1,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_floor), numphases);
                z_y(2,1,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_ceil), numphases);
                z_y(2,2,1) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_floor), numphases);
                z_y(2,2,2) = v_y(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_ceil), numphases);
            
                z_z(1,1,1) = v_z(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_floor), numphases);
                z_z(1,1,2) = v_z(uint16(y_temp_floor), uint16(x_temp_floor), uint16(z_temp_ceil), numphases);
                z_z(1,2,1) = v_z(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_floor), numphases);
                z_z(1,2,2) = v_z(uint16(y_temp_floor), uint16(x_temp_ceil), uint16(z_temp_ceil), numphases);
                z_z(2,1,1) = v_z(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_floor), numphases);
                z_z(2,1,2) = v_z(uint16(y_temp_ceil), uint16(x_temp_floor), uint16(z_temp_ceil), numphases);
                z_z(2,2,1) = v_z(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_floor), numphases);
                z_z(2,2,2) = v_z(uint16(y_temp_ceil), uint16(x_temp_ceil), uint16(z_temp_ceil), numphases);
                %-------------------------------------------------------
                v_x_temp = interp3([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil],[z_temp_floor z_temp_ceil], z_x, px, py, pz,'linear');
                v_y_temp = interp3([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil],[z_temp_floor z_temp_ceil], z_y, px, py, pz,'linear');
                v_z_temp = interp3([x_temp_floor x_temp_ceil], [y_temp_floor y_temp_ceil],[z_temp_floor z_temp_ceil], z_z, px, py, pz,'linear');
                
                if isnan(v_x_temp)
                    v_x_temp=0;
                end
                if isnan (v_y_temp)
                    v_y_temp=0;
                end
                if isnan (v_z_temp)
                    v_z_temp=0;
                end
                
                %-------------------------------------------------------
                vx(j,numphases) = v_x_temp;
                vy(j,numphases) = v_y_temp;
                vz(j,numphases) = v_z_temp;
                vr(j,numphases) = sqrt(v_x_temp^2 + v_y_temp^2+v_z_temp^2);
                
            
                
            end % last phase velocity calculation
        end % j loop (# of seeds)
    end % if
end % i loop (# of phases)

multiWaitbar('Tracking...', 'Close');

