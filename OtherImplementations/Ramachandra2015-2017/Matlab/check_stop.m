function [flg_ts, flg_n] = check_stop(t, t_stop, dt)

if t>t_stop
    flg_ts = 0;
else
    flg_ts = 1;
end

flg_n = 1;
