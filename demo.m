%Demo - Run this code to have a demonstration how to use MOESP -tool. By
%applying some random simulation data, and LTI representation of
%input-output behaviour is found

%NOTE: The data file containing the random data must be stored in same
%folder for this demo to run.
clear;clc;

%% Load the IO-data to estimate:
load('someData.mat','inputs','outputs');
%Note: the inputs and outputs must be in "deviation from nominal" format.

%% Apply the MOESP identification:
[A,B,C,D,~,~] = MOESP(inputs,outputs);

%% Visually evaluate how the model fits the original data:
x = zeros(size(A,1),1); %Initial state is zero deviation from nominal.
y = zeros(size(outputs)); %Initialize data saving
for k = 1:size(inputs,2) %Simulate through the data
    x = A*x + B*inputs(:,k); %No need to save state, so overwrite
    y(:,k) = C*x + D*inputs(:,k); %Save the obtained LTI estimate output
end
for iOutput = 1:size(outputs,1)
    figure();plot(outputs(iOutput,:));
    hold on;plot(y(iOutput,:)); hold off
    legend('Original measurement','LTI estimate')
end