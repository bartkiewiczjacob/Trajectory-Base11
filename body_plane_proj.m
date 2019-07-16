function projection = body_plane_proj(vector, ref)
% projects vector on the plane perpendicular to ref

if ref == 'y'
    projection = vector - dot(vector,[0; 1; 0])*[0; 1; 0];
elseif ref == 'z'
    projection = vector - dot(vector,[0; 0; 1])*[0; 0; 1];
else
    error('Enter 'y' or 'z');
end

end