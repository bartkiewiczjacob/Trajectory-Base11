% Andrew Meyer
% 06/16/2019

% Quaternion Multiplication Function

function p = Quat_Mult(q1, q2)
    q1_s = q1(1);
    q1_v = q1(2:4);
    q2_s = q2(1);
    q2_v = q2(2:4);
    p_s = q1_s*q2_s - dot(q1_v,q2_v);
    p_v = q1_s*q2_v + q2_s*q1_v + cross(q1_v,q2_v);
    [rows,~] = size(p_v);
    if rows == 1
        p = [p_s, p_v];
    else
        p = [p_s; p_v];
    end
end