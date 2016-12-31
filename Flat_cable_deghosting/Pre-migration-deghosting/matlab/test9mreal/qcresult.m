close all
clear all
clc

trid = 280;

fid = fopen('filter9mghost.bin','r');
datain = fread(fid,[4000 576],'single');
fclose(fid);
fid = fopen('output_01_50_5s.bin','r');
dataout = fread(fid,[4000 576],'single');
fclose(fid);
t = [0:0.002:0.002*3999]';
plot(t,datain(:,trid),'red',t,dataout(:,trid),'blue');
