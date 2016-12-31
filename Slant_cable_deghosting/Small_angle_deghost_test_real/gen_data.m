clear all 
close all
clc

nt=2001;
nx=648;

fid = fopen('Tam_data.bin','rb');
input = fread(fid,[nt nx],'float');
fclose(fid);


smallinput = input(:,1:300);

fid = fopen('Small_input.bin','wb');
fwrite(fid, smallinput,'single');
fclose(fid);


