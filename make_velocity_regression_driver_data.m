%% make velocity regression data 
function [data_use] = make_velocity_regression_driver_data(driver_data)


speed = driver_data(:,5);
speed_shifted = [0; speed];
speed_shifted_cut = speed_shifted(2: size(speed_shifted,1)-1);
speed_use = speed_shifted_cut;
acc = driver_data(:,6);
acc_use = acc(2:size(acc,1));
speed_frontcar = driver_data(:,9);
speed_frontcar_use = speed_frontcar(2:size(speed_frontcar,1));
acc_frontcar= driver_data(:,10);
acc_frontcar_use = acc_frontcar(2:size(acc_frontcar,1));
range = driver_data(:,7);
range_use = range(2:size(range,1));
range_rate = driver_data(:,8);
range_rate_use = range_rate(2:size(range_rate,1));
kdb = driver_data(:,11);
kdb_use = kdb(2:size(kdb,1));
jerk = driver_data(:,12);
jerk_use = jerk(2:size(jerk,1));
invTTC = driver_data(:,13);
invTTC_use = invTTC(2:size(invTTC,1));
THW = driver_data(:,14);
THW_use = THW(2:size(THW,1));

data_use = [speed_use(:,1) acc_use(:,1) speed_frontcar_use(:,1) acc_frontcar_use(:,1) range_use(:,1) range_rate_use(:,1) kdb_use(:,1) jerk_use(:,1) invTTC_use(:,1) THW_use(:,1)];

end
