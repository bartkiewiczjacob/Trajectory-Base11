function projection = body_plane_proj(vector, ref)
% projects vector on the plane perpendicular to ref

if ref == 'y'
    projection = [vector(1); 0; vector(3)];
elseif ref == 'z'
    projection = [vector(1); vector(2); 0];
else
    error('Enter "y" or "z"');
end

end