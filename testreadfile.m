fid = fopen('../MyVTKData/Pigrgba.bin','r','l');
rgba = fread(fid, [4,512*512*134], 'float','l')';
%rgba=reshape(rgba,512,512,134);
fclose(fid);