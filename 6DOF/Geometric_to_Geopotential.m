% Andrew Meyer
% 06/18/2019

function geopotential_alt = Geometric_to_Geopotential(geometric_alt)
    R = 6371000; % m [1]
    geopotential_alt = (R / (R + geometric_alt)) ...
                                   * geometric_alt; % m [2] 
end

% References
% [1] Earth's volumetric mean radius:
%     https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
% [2] Relation between geopotential and geometric altitude
%     Anderson, J. "Introduction to Flight", 8th Edition, pg. 116