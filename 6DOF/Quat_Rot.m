% Andrew Meyer
% 06/15/2019

% Quaternion Rotation Function

function r_rot = Quat_Rot(r, q)
    q1 = q/norm(q);
    q1_s = q1(1);
    q1_v = q1(2:4);
    q2_s = q1_s;
    q2_v = -q1_v;
    p1_s = -dot(q1_v,r);
    p1_v = q1_s*r + cross(q1_v,r);
    %p2_s = p1_s*q2_s - dot(p1_v,q2_v);
    %p2_v = p1_s*q2_v + q2_s*p1_v + cross(p1_v,q2_v);
    r_rot = p1_s*q2_v + q2_s*p1_v + cross(p1_v,q2_v);
end