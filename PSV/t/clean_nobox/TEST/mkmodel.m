vp=zeros(1024,1024)+6;
vs=zeros(1024,1024)+4;
den=zeros(1024,1024)+3.0;

fid=fopen('homo.vp','wb');
fwrite(fid,vp,'single');
fclose(fid);

fid=fopen('homo.vs','wb');
fwrite(fid,vs,'single');
fclose(fid);

fid=fopen('homo.den','wb');
fwrite(fid,den,'single');
fclose(fid);
