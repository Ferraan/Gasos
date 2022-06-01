clear;
data=csvread('Treball_sortida.csv');
data(1,:)=[];
%plot(data(:,2),data(:,3));
%plot(data(:,2),data(:,4));
plot(data(1:end-1,2),data(1:end-1,13));
