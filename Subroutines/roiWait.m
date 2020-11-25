function mask = roiWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@roiClickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
mask = createMask(hROI);

end