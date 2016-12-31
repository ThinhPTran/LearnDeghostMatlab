% Test leastsub

clear all
close all
clc

tic

epsilon = 100.0;
iter_cg = 1000;

nt = 500;
nx = 401;


fid = fopen('upcon_ghost1','rb');
down = fread(fid,[nt nx],'single');
fclose(fid);

fid = fopen('receiver_data_p','rb');
primary = fread(fid,[nt nx],'single');
fclose(fid);

error = zeros(nt,nx);


error = leastsub(primary,down, nx, nt, epsilon, iter_cg);


fid = fopen('error','wb');
fwrite(fid, error,'single');
fclose(fid);


toc