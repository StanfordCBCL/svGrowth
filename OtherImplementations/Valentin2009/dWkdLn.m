function  y = dWkdLn(Ln, c_c2, c_c3, c_c2_c, c_c3_c)

if Ln >= 1
    y = c_c2*Ln*(Ln^2 - 1)*exp(c_c3*(Ln^2 - 1)^2);
else
    y = c_c2_c*Ln*(Ln^2 - 1)*exp(c_c3_c*(Ln^2 - 1)^2);
end