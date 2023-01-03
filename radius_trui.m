function   r = radius_trui(xy) % Return grading grain radius at location (x,y)
persistent F

if isempty(F)                           % Read in the trui image once
    A  = double(imread('trui.png','PNG'));
    A = flipud(A(:,:,1)); 
    rf = @(s) 0.002+0.006*s+0.012*s.^8; % Conversion of brightness to grain radius 
    F  = rf(A/255);
end

ixy = round(255*xy);                      % Given a location (x,y), evaluate  
r = F(ixy(:,2)+256*ixy(:,1)+1);           % the corresponding grain radius