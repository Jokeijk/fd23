function [b,h]=isisread(file)
%file='W.isis';
at=dir(file);
fid=fopen(file,'rb');
%fid=fopen('/home/clay/FD/acou/stag/FaultBendFold.isis','r');
a=fread(fid,8,'int32');
nt=a(8);
nrec=at.bytes/(nt+21)/4;
fclose(fid);
fid=fopen(file,'rb');
a=fread(fid,(nt+21)*nrec,'float=>float');
a=reshape(a,nt+21,nrec);
b=a(22:nt+21,:);
if nargout==2
    h.xs=a(1,:)';
    h.zs=a(3,:)';
    h.xr=a(4,:)';
    h.zr=a(6,:)';
    h.dt=a(9,:)';
    h.gain=a(13,:)';
end
fclose(fid);
