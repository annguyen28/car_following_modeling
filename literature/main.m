%% this codes simulates car following with developed model and gipps model

close all
clear all

%% load needed variables
load('a_n.mat')
load('b_n.mat')
load('lead_b.mat')
load('S.mat')
load('T.mat')
load('V_n.mat')

load('max_speed_a.mat')
load('max_speed_b.mat')
load('max_range.mat')
load('max_THW.mat')
load('max_TTC_inverse.mat')
load('max_range_r.mat')
load('max_lead_acc.mat')
load('max_jerk')

load ('para_Sim_mode_1.mat')
load ('para_Sim_mode_2.mat')
load ('para_Sim_mode_3.mat')
load ('para_Sim_mode_4.mat')

load('pattern1.mat')
load('pattern2.mat')
load('pattern3.mat')
load('pattern4.mat')
load('pattern5.mat')
load('pattern6.mat')
load('pattern7.mat')

%% adding points to leading car velocity pattern
% pattern3AddTime = (max(pattern3(:,1))+0.1):0.1:((max(pattern3(:,1))+(length(pattern3)*0.1)));
% pattern3AddMore = [pattern3AddTime.',flip(pattern3(:,2:end))];
% pattern = [pattern3;pattern3AddMore];
pattern = pattern2;
t = 0.3;
simulationPeriod = 0.3;
modelResponseTime = t/simulationPeriod;
time = (min(pattern(:,1)):simulationPeriod:max(pattern(:,1))).';
lead_velo = pchip(pattern(:,1),pattern(:,2),time);
lead_range = pchip(pattern(:,1),pattern(:,3),time);
lead_thw = lead_range(:)./lead_velo(:);
time = time - min(time);
initialDistance = 10; %[m]

%% Initialization of the variables to stack the results
position = zeros(length(time), 3);
velocity = zeros(length(time), 3);
acceleration = zeros(length(time), 3);

range     = zeros(length(time), 2);
rangeRate = zeros(length(time), 2);
jerk = zeros(length(time), 2);
THW = zeros(length(time), 1);
TTC_inverse = zeros(length(time), 1);
modelOutput = zeros(length(time), 2);
modeClassifier = [11, 21, 31];

%% Main Loop
% lead car initialization
position(1,1) = initialDistance;
velocity(:,1) = lead_velo(:,1);
acceleration(1,1) = 0;

% ego car our model and gipps model
position(1,2:3) = position(1,1) - initialDistance;
velocity(1,2:3) = lead_velo(1,1);
acceleration(1,2:3) = 0;
jerk(1,:) = 0;

modelOutput(1,:) = velocity(1,2:3);

range(1,:) = initialDistance;
rangeRate(1,:) = velocity(1,1) - velocity(1,2:3);
THW(1) = range(1,1)/velocity(1,2);
TTC_inverse(1) = (-rangeRate(1,1)/range(1,1));

%%

for i = 2:length(time)
    
    % update leading car state
    position(i,1)      = position(i-1,1) + (velocity(i-1,1)*simulationPeriod);
    acceleration(i,1)  = (velocity(i,1) - velocity(i-1,1))/simulationPeriod;
    
    % ego car: our model and gipps model
    if time(i) <= t
        position(i,2:3) = position(i-1,2:3) + velocity(i-1,2:3)*simulationPeriod;
        velocity(i,2:3) = lead_velo(i,1);
        acceleration(i,2:3) = (velocity(i,2:3) - velocity(i-1,2:3))/simulationPeriod;
        range(i,:) = position(i,1) - position(i,2:3);
        rangeRate(i,:) = velocity(i,1) - velocity(i,2:3);
        jerk(i,:) = (acceleration(i,2:3) - acceleration(i-1,2:3))/simulationPeriod;
        THW(i) = range(i,1)/velocity(i,2);
        TTC_inverse(i)= (-rangeRate(i,1))/range(i,1);
        modelOutput(i,:) = velocity(i,2:3);
    else
        % gipps model
        position(i,3) = position(i-1,3) + velocity(i-1,3)*simulationPeriod;
        gippsTmp1 = velocity(i-modelResponseTime,3) + (2.5 * a_n * T * ...
            (1 - (velocity(i-modelResponseTime,3)/V_n)) * ...
            (sqrt(0.025 + (velocity(i-modelResponseTime,3)/V_n))));
        gippsTmp2 = b_n * T + (sqrt(abs((b_n^2*T^2) - ...
            b_n * (2*(range(i-modelResponseTime,2)-S)-...
            (velocity(i-modelResponseTime,3)*T)-...
            (velocity(i-modelResponseTime,1)^2/lead_b)))));
        velocity(i,3) = min(gippsTmp1, gippsTmp2);
        acceleration(i,3)  = (velocity(i,3) - velocity(i-1,3))/simulationPeriod;
        range(i,2) = position(i,1) - position(i,3);
        rangeRate(i,2) = velocity(i,1) - velocity(i,3);
        jerk(i,2) = (acceleration(i,3) - acceleration(i-1,3))/simulationPeriod;
        modelOutput(i,2) = velocity(i,3);
        
        % our model
        position(i,2) = position(i-1,2) + velocity(i-1,2)*simulationPeriod;
        temp_x = [1, velocity(i-modelResponseTime,2)/max_speed_b, range(i-modelResponseTime,1)/max_range,...
                  rangeRate(i-modelResponseTime,1)/max_range_r, jerk(i-modelResponseTime,1)/max_jerk,...
                  TTC_inverse(i-modelResponseTime)/max_TTC_inverse, THW(i-modelResponseTime)/max_THW,...
                  acceleration(i-modelResponseTime,1)/max_lead_acc];
        if(range(i,1) < modeClassifier(1))
            velocity(i,2) = temp_x * para_Sim_mode_1 * max_speed_a;
        elseif(range(i-1,1)>= modeClassifier(1) && range(i-1,1) < modeClassifier(2))
            velocity(i,2) = temp_x * para_Sim_mode_2 * max_speed_a;
        elseif(range(i-1,1)>= modeClassifier(2) && range(i-1,1) < modeClassifier(3))
            velocity(i,2) = temp_x * para_Sim_mode_3 * max_speed_a;
        else
            velocity(i,2) = temp_x * para_Sim_mode_4 * max_speed_a;
        end
        acceleration(i,2)  = (velocity(i,2) - velocity(i-1,2))/simulationPeriod;
        range(i,1) = position(i,1) - position(i,2);
        rangeRate(i,1) = velocity(i,1) - velocity(i,2);
        jerk(i,1) = (acceleration(i,2) - acceleration(i-1,2))/simulationPeriod;
        THW(i) = range(i,1)/velocity(i,2);
        TTC_inverse(i)= (-rangeRate(i,1))/range(i,1);
        
        modelOutput(i,1) = velocity(i,2);
        
    end
    
end

model_thw = THW;
Gipps_thw = range(:,2)./velocity(:,3);

thw = [model_thw, Gipps_thw];


%% plotting results 
figure(1)
hold on 
grid on 
grid minor 
plot(velocity*3.6)
legend('Leading car', 'proposed model', 'Gipps model')
xlabel('time [s]')
ylabel('speed [km/h]')
title('Car following simulation for driver 4 as leading car, speed profiles')
hold off 

figure(2)
hold on 
grid on 
grid minor 
plot(acceleration)
legend('Leading car', 'proposed model', 'Gipps model')
xlabel('time [s]')
ylabel('acceleration [m/s^2]')
title('Car following simulation for driver 4 as leading car, acceleration profiles')
hold off 

figure(3)
hold on 
grid on 
grid minor 
plot(range)
legend('proposed model', 'Gipps model')
xlabel('time [s]')
ylabel('range [m]')
title('Car following simulation for driver 4 as leading car, distance to lead car')
hold off 

figure(4)
hold on 
grid on 
grid minor 
plot(thw)
legend('proposed model', 'Gipps model')
xlabel('time [s]')
ylabel('Time headway [s]')
title('Car following simulation for driver 4 as leading car, Time headway')
hold off 