function y = f_beta(beta, kq)

% gamma = 0.1;
% 
% y_f = 5.0; % was 5.0
% if beta < y_f
%     y = 0.0;
% else
% %	y = (beta - y_f)^2 * kq ;
% % y = (1 + abs(beta - y_f)) * kq;
% 	y = 9*kq*(1 - exp(-gamma * (beta - y_f)^2));
% end

y = abs(beta - 1) * kq;