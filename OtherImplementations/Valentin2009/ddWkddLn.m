function  y = ddWkddLn(Ln, c_c2, c_c3, c_c2_c, c_c3_c)

if Ln >= 1
    y = c_c2*(3*Ln^2 - 1 + 4*c_c3*(Ln*(Ln^2 - 1))^2)*exp(c_c3*(Ln^2 - 1)^2);
else
    exp_Q = exp(c_c3_c*(Ln^2-1.0)^2);
    y = c_c2_c*(3*Ln^2-1.0+4*c_c3_c*(Ln*(Ln^2-1.0))^2)*exp_Q;
end