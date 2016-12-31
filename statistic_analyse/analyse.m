clear all
close all
clc

data = load('statistic_z2.txt','r'); 


size(data)

src_depth = data(:,2); 
rec_depth = data(:,3); 


figure(); 
hist(src_depth); 


figure(); 
hist(rec_depth); 









