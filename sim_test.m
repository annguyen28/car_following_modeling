%% load filtered data
load ("driver1_data1.mat");
load ("driver1_data2.mat");
load ("driver2_data1.mat");
load ("driver2_data2.mat");
load ("driver3_data1.mat");
load ("driver3_data2.mat");
load ("driver4_data1.mat");
load ("driver4_data2.mat");
load ("driver5_data1.mat");
load ("driver5_data2.mat");
load ("driver6_data1.mat");
load ("driver6_data2.mat");
load ("driver7_data1.mat");
load ("driver7_data2.mat");

input_data = [driver1_data1(:,:)];
%     ;driver2_data1(:,:);driver3_data1(:,:)];
%     ;driver4_data1(:,:);driver5_data1(:,:);driver6_data1(:,:);driver7_data1(:,:)];
% driver1_data2(:,:);driver2_data2(:,:);driver3_data2(:,:);;driver4_data2(:,:);driver5_data2(:,:);driver6_data2(:,:);driver7_data2(:,:)];
driver_data = make_velocity_regression_driver_data(input_data); 
% driver_data = [driver1_data1(:,:)];           
%% variables in data

% '1   time';
% '2   throttle';
% '3   brake'; 
% '4   steer'; 
% '5   speed'; 
% '6   acceleration'; 
% '7   range'; 
% '8   range_rate';
% '9   lead_car_velocity';
% '10  lead_car_acc';
% '11  kdb';
% '12  jerk';
% '13  TTC_inverse';
% '14  THW' ];

%% nomalize data

% n_time = driver_data(:,1)/max(abs(driver_data(:,1)));
% n_throttle = driver_data(:,2)/max(abs(driver_data(:,2)));
n_speed = 2*((driver_data(:,1)-min(driver_data(:,1)))/(max(driver_data(:,1))-min(driver_data(:,1))))-1;
n_acceleration = 2*((driver_data(:,2)-min(driver_data(:,2)))/(max(driver_data(:,2))-min(driver_data(:,2))))-1;
n_range = 2*((driver_data(:,5)-min(driver_data(:,5)))/(max(driver_data(:,5))-min(driver_data(:,5))))-1;
n_range_rate = 2*((driver_data(:,6)-min(driver_data(:,6)))/(max(driver_data(:,6))-min(driver_data(:,6))))-1;
n_lead_car_velocity = 2*((driver_data(:,3)-min(driver_data(:,3)))/(max(driver_data(:,3))-min(driver_data(:,3))))-1;
n_lead_car_acc = 2*((driver_data(:,4)-min(driver_data(:,4)))/(max(driver_data(:,4))-min(driver_data(:,4))))-1;
n_kdb = 2*((driver_data(:,7)-min(driver_data(:,7)))/(max(driver_data(:,7))-min(driver_data(:,7))))-1;
n_jerk = 2*((driver_data(:,8)-min(driver_data(:,8)))/(max(driver_data(:,8))-min(driver_data(:,8))))-1;
n_TTC_inverse = 2*((driver_data(:,9)-min(driver_data(:,9)))/(max(driver_data(:,9))-min(driver_data(:,9))))-1;
n_THW = 2*((driver_data(:,10)-min(driver_data(:,10)))/(max(driver_data(:,10))-min(driver_data(:,10))))-1;


%% assign y - phi
y  = n_acceleration;
% offset_matrix = ones(length(y),1);
phi = [n_speed n_range n_range_rate n_lead_car_velocity n_lead_car_acc n_kdb n_jerk n_TTC_inverse n_THW];

%%
clear driver1_data1 driver2_data1 driver3_data1 driver4_data1 driver5_data1 driver6_data1 driver7_data1
clear driver1_data2 driver2_data2 driver3_data2 driver4_data2 driver5_data2 driver6_data2 driver7_data2
%% estimate feature vectors

opt_f.c = 1000;
opt_f.rmv_const = true; 

% following four options are default settings of calculation.
% opt_f.calc_r = true;  
% opt_f.calc_ir = true;
% opt_f.calc_spr = true;
% opt_f.calc_w = true;

% calculate the feature vectors through the dynamics.
[gLDs, LDs] = ohpk_pwarx_data2feature_space( phi, y, opt_f );
figure;plotmatrix(gLDs);
%%
figure
E = evalclusters([gLDs(:,1),gLDs(:,11)],'kmeans','DaviesBouldin','KList',[1:10]);
plot(E)

%% clustering

mode_num = 2;
%mode_num = E.OptimalK
opt.NumOfInitialValues = 500;   
% opt.MaxRepetations = 500;       
opt.CenterInitializeMethod = 'pickout';    
% opt.CenterInitializeMethod = 'normal';    
% opt.CenterInitializeMethod = 'uniform';    
% opt.CenterInitializeStd = std(gLDs);    
% opt.CenterInitializeMean = mean(gLDs,1);    

tic;
opt.ShowProgress = 't';
opt.ShowProgressSkip = 100;
[center, class] = ohpk_pwarx_weighted_kmeans(gLDs, mode_num, LDs, opt);
toc;

%%
%perate mode
for i = 1:mode_num
    data(i).ymode = y(class==i); 
    data(i).phimode = phi(class==i, :); 
end

%%
 %mode 1

% create binary table
number_of_variable = size(data(1).phimode,2);
binary_table = dec2bin(1:(2^(number_of_variable)-1)) - '0';

for i = 1:1:2^(number_of_variable)-1;
data_table1(i).mode_index = i;
data_table1(i).binary_table = [binary_table(i,:)];
end
 [row col] = find(binary_table);
  %calculate aic
for i= 1:1:size(binary_table,1);
   selected_array1 = data(1).phimode(:,col(row == i));
   data_table1(i).model_variables_matrices = selected_array1(:,:);
   [data_table1(i).theta,~,~,~,data_table1(i).logL] = mvregress(data_table1(i).model_variables_matrices ,data(1).ymode);
   aicmode1(i) = aicbic(data_table1(i).logL, size(data_table1(i).model_variables_matrices,2));
   invaicmode1(i) = 1/aicmode1(i);
end
[aic_value_mode1 estimated_model_index_mode1] = min(aicmode1);

beta1 = inv(data(1).phimode'*diag(LDs.w(class==1))*data(1).phimode)*data(1).phimode'*diag(LDs.w(class==1))*data(1).ymode;

fprintf('\nVaribles of mode 1 are: \n');
fprintf(' %d  ',data_table1(estimated_model_index_mode1).binary_table);
fprintf('\n');
fprintf('With parameters: \n');
%fprintf(' %f  ',data_table1(estimated_model_index_mode1).theta);
fprintf(' %f  ',beta1);
fprintf('\n');
%%
%mode 2

% create binary table
number_of_variable = size(data(2).phimode,2);
binary_table = dec2bin(1:(2^(number_of_variable)-1)) - '0';

for i = 1:1:2^(number_of_variable)-1;
data_table2(i).mode_index = i;
data_table2(i).binary_table = [binary_table(i,:)];
end
 [row col] = find(binary_table);
   

%calculate aic
for i= 1:1:size(binary_table,1);
   selected_array2 = data(2).phimode(:,col(row == i));
   data_table2(i).model_variables_matrices = selected_array2(:,:);
   [data_table2(i).theta,~,~,~,data_table2(i).logL] = mvregress(data_table2(i).model_variables_matrices ,data(2).ymode);
   aicmode2(i) = aicbic(data_table2(i).logL, size(data_table2(i).model_variables_matrices,2));
end
[aic_value_mode2 estimated_model_index_mode2] = min(aicmode2);

beta3 = pinv(data(2).phimode'*diag(LDs.w(class==2))*data(2).phimode)*data(2).phimode'*diag(LDs.w(class==2))*data(2).ymode;

fprintf('\nVaribles of mode 3 are: \n');
fprintf(' %d  ',data_table2(estimated_model_index_mode2).binary_table);
fprintf('\n');
fprintf('With parameters: \n');
%fprintf(' %f  ',data_table2(estimated_model_index_mode2).theta);
fprintf(' %f  ',beta3);
fprintf('\n');

%%
 %mode 3

% create binary table
number_of_variable = size(data(3).phimode,2);
binary_table = dec2bin(1:(2^(number_of_variable)-1)) - '0';

for i = 1:1:2^(number_of_variable)-1;
data_table3(i).mode_index = i;
data_table3(i).binary_table = [binary_table(i,:)];
end
 [row col] = find(binary_table);
   

%calculate aic
for i= 1:1:size(binary_table,1);
   selected_array3 = data(3).phimode(:,col(row == i));
   data_table3(i).model_variables_matrices = selected_array3(:,:);
   [data_table3(i).theta,~,~,~,data_table3(i).logL] = mvregress(data_table3(i).model_variables_matrices ,data(3).ymode);
   aicmode3(i) = aicbic(data_table3(i).logL, size(data_table3(i).model_variables_matrices,2));
end
[aic_value_mode3 estimated_model_index_mode3] = min(aicmode3);

beta2 = pinv(data(3).phimode'*diag(LDs.w(class==3))*data(3).phimode)*data(3).phimode'*diag(LDs.w(class==3))*data(3).ymode;
fprintf('\n Varibles of mode 2 are: \n');

fprintf(' %d  ',data_table3(estimated_model_index_mode3).binary_table);
fprintf('\n');
fprintf('With parameters: \n');
%fprintf(' %f ',data_table3(estimated_model_index_mode3).theta);
fprintf(' %f  ',beta2);
fprintf('\n');



%%
% plot data after cluster
clr = lines(mode_num);
figure('name', 'Clustered Data'), hold on
% scatter3(phi(:,1)', phi(:,2)',y(:,:)', 10, 'Marker','d','LineWidth',1);
scatter3(data(2).phimode(:,1), data(2).phimode(:,4),data(2).ymode(:,:), 'Marker','d', 'LineWidth',1)
hold off
view(3), axis vis3d, box on, rotate3d on
xlabel('x(1)'), ylabel('x(2)'), zlabel('y');
grid on
set(gca,'Fontsize',11)
set(gcf,'Position',[10 10 600 400])
%%
% plot data after cluster
clr = lines(mode_num);
figure('name', 'Clustered Data'), hold on
scatter3(phi(:,2), phi(:,3),y(:,1),10, clr(class,:), 'Marker','x','LineWidth',1);
% scatter3(gLDs(:,6), gLDs(:,15),y(:,1),50, clr(class,:), 'Marker','.','LineWidth',1)
hold off
view(3), axis vis3d, box on, rotate3d on
% xlabel('speed'), ylabel('kdb'), zlabel('acc');

view([0 0]);grid on
set(gca,'Fontsize',14)
set(gcf,'Position',[10 10 600 400])

%%
% plot data after cluster
clr = lines(mode_num);
figure('name', 'Clustered Data'), hold on
scatter3(gLDs(:,3), gLDs(:,18),gLDs(:,6),50, clr(class,:), 'Marker','.','LineWidth',1)
% scatter3(center(:,1), center(:,2), center(:,3), 100, clr, 'Marker','o', 'LineWidth',1)
hold off
view([0 0]), axis vis3d, box on, rotate3d on
%xlabel('range'), ylabel('kdb'), zlabel('jerk');
grid on
% set(gca,'Fontsize',13)
% set(gcf,'Position',[50 50 600 400])
% clear gLDs n_thottle



%%
cost_func = "NRMSE";
yref = n_speed(class==1);
beta = pinv(data(2).phimode'*data(2).phimode)*data(2).phimode'*data(2).ymode
Y = data(2).phimode*[0 0 -1.2 -7 0 3]';
Ynew = data(2).phimode*beta3;
fit1 = goodnessOfFit(y3, y3ref,cost_func)
fit2 = goodnessOfFit(y3new, y3ref,cost_func)



plot([1:1:size(yref)],y3ref,'Marker','o', 'LineWidth',2,'Color','r')
hold on
plot([1:1:size(yref)],y3,'Marker','diamond', 'LineWidth',1,'Color','y')
plot([1:1:size(yref)],y3new,'Marker','x', 'LineWidth',0.1, 'Color','k')
