clear all
close all
clc

nt = 500;
nx = 401;


fid = fopen('upcon_ghost1','rb');
down = fread(fid,[nt nx],'single');
fclose(fid);

fid = fopen('receiver_data_p','rb');
primary = fread(fid,[nt nx],'single');
fclose(fid);


err = abs(primary-down);
fid = fopen('error','wb');
fwrite(fid, err,'single');
fclose(fid);




