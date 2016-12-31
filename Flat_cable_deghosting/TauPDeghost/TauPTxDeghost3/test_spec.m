clear all 
close all
clc


nt = 500;
nx = 401;


fid = fopen('data_with_ghost','rb');
input = fread(fid,[nt nx],'single');
fclose(fid);

input_fk = abs(fft2(input));


fid = fopen('deghost.bin','rb');
output = fread(fid,[nt nx],'single');
fclose(fid);

output_fk = abs(fft2(output));


figure();
imagesc(input_fk);


figure();
imagesc(output_fk);


fid = fopen('input_fk','wb');
fwrite(fid,input_fk,'single');
fclose(fid);


fid = fopen('output_fk','wb');
fwrite(fid,output_fk,'single');
fclose(fid);





