close all
clear all
clc

dt= 0.00066666
t = [0:dt:3000*dt]';


fid = fopen('flat_primary_160.bin','r');
primary160 = fread(fid,[3001 100],'single');
fclose(fid);

fid = fopen('flat_primary_100.bin','r');
primary100 = fread(fid,[3001 100],'single');
fclose(fid);

fid = fopen('flat_ghost.bin','r');
prighost = fread(fid,[3001 100],'single');
fclose(fid);

fid = fopen('flat_ghost_out.bin','r');
prighostout = fread(fid,[3001,100],'single');
fclose(fid);


figure(1)
plot(t,primary160(1:3001,40),'green',t,primary100(1:3001,40),'yellow',t,prighost(1:3001,40),'blue',t,prighostout(1:3001,40),'red');
legend('primary160','primary100','ghost data','deghosted data');