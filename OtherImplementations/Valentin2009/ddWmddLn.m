function  y = ddWmddLn(Ln, c_m2, c_m3, c_m2_c, c_m3_c)

if Ln >= 1
    y = c_m2*(3*Ln^2 - 1 + 4*c_m3*(Ln*(Ln^2 - 1))^2)*exp (c_m3*(Ln^2 - 1)^2);
else
    exp_Q = exp (c_m3_c*(Ln^2-1.0)^2);
    y=c_m2_c*(3*Ln^2-1.0+4*c_m3_c*(Ln*(Ln^2-1.0))^2)*exp_Q;
end
