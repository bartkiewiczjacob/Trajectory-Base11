% AUTHOR: Andy Meyer
% LAST MODIFIED: 11/07/2019

% FUNCTION PURPOSE:
% Compute the inverse of a quaternion.

% INPUT:
% q

% OUTPUT:
% q^-1

function q_inverse = Quat_Inv(q)
    q_s = q(1);
    q_v = q(2:4);
    q_norm = norm(q);
    q_inverse = [q_s, -q_v] ./ q_norm;
end