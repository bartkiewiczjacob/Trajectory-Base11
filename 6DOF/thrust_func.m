% AUTHOR: Andy Meyer
% LAST MODIFIED: 11/07/2019

% FUNCTION PURPOSE
% Compute the thrust of the rocket based on the time from t-zero

% INPUTS
% t - time in seconds from t-zero

% OUTPUTS
% thrust - thrust in Newtons

function thrust = thrust_func(t)
    t_burnout = 9; % s
    if t < t_burnout
        thrust = 900*4.44822; % N
    else
        thrust = 0; % N
    end
end