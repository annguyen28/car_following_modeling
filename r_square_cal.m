%% r square

y_mode1= data(1).ymode;
sum_var_square=sum((y_mode1-mean(y_mode1)).^2);
e_square_sum = sum((data(1).ymode- data(1).phimode*beta1).^2);
r_square = 1 - e_square_sum/sum_var_square

cost_func = "NRMSE";
y_ref = data(1).ymode;
% beta = pinv(data(2).phimode'*data(2).phimode)*data(2).phimode'*data(2).ymode
y_est = data(1).phimode*([0.00722301331497579;0.0359878428817944;0.0217964776228112;0.0625458271114140;0.177143217264628;-0.00427534179047799;0.0415963176839847;-1.11054224540449;0.0983182082264551]);

fit1 = goodnessOfFit(y_ref, y_est,cost_func)

