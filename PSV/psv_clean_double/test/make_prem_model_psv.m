%% make model for fd
prem=mkreadclr([mkttboxdata 'prem.clr'],'silent');

[~,moho_flat]=mksfer2flat(-1,prem.moho,prem.rp);
[~,bottom_flat]=mksfer2flat(-1,5200,prem.rp);

h=moho_flat/15.5;

nz=round((bottom_flat/h)/32+3)*32;

z_f=(0:nz-1)*h; 
[~,z_s]=mkflat2sfer(-1,z_f,prem.rp);

z_f=z_f(:);
z_s=z_s(:);

vp_s=z_s*0;
vs_s=z_s*0;
rho_s=z_s*0;

for i=1:length(z_s)
    [~,vp_s(i)]=mksampleclr(prem,z_s(i),'vp');
    [~,vs_s(i)]=mksampleclr(prem,z_s(i),'vs');
    [~,rho_s(i)]=mksampleclr(prem,z_s(i),'rho');
end
%% no crust
id=find(z_s>24.4,1)-1;
vp_s(1:id)=mean(vp_s(1:id));
vs_s(1:id)=mean(vs_s(1:id));
rho_s(1:id)=mean(rho_s(1:id));
%%


[vp_f,~]=mksfer2flat(vp_s,z_s,prem.rp);
[vs_f,~]=mksfer2flat(vs_s,z_s,prem.rp);

% m=3 is exact for SH
% m=0 may be more accurate for Reighleigh
% m=-2
m=-2;

rho_f=rho_s.*(prem.rp./(prem.rp-z_s)).^-(m+2);
%%
%% New Depth
evdp=39.3;
[~,ievdp]=min(abs(evdp-z_s));
evdp_fd=z_s(ievdp);
evdp_flat=z_f(ievdp);
fid=fopen('source_pos','w');
fprintf(fid,'evdp=%f\n',evdp);
fprintf(fid,'evdp_fd=%f\n',evdp_fd);
fprintf(fid,'evdp_flat=%f\n',evdp_flat);
fprintf(fid,'ievdp=%d\n',ievdp);
fprintf(fid,'vs_f=%f\n',vs_f(ievdp));
fprintf(fid,'vp_f=%f\n',vp_f(ievdp));
fprintf(fid,'rho_f=%f\n',rho_f(ievdp));
fclose(fid);
%% FK model
fid=fopen('prem_psv_fk','w');
tmp=[vs_f*0+h,vs_f, vp_f, rho_f, vs_f*0+1E10, vs_f*0+1E10];
tmp(end,1)=0;
tmp(1,1)=h/2;

fprintf(fid,'%f %f %f %f %f %f\n',tmp');
fclose(fid);
%% FD model

nx=round(deg2km(150)/h/32)*32;

toppad=640;
vs_f=[vs_f(1)+zeros(toppad,1);vs_f];
tmp=repmat(vs_f,1,nx); tmp=tmp';
fid=fopen('prem_psv.vs','wb');
fwrite(fid,tmp,'single');
fclose(fid);

vp_f=[vp_f(1)+zeros(toppad,1);vp_f];
tmp=repmat(vp_f,1,nx); tmp=tmp';
fid=fopen('prem_psv.vp','wb');
fwrite(fid,tmp,'single');
fclose(fid);

rho_f=[rho_f(1)+zeros(toppad,1);rho_f];
tmp=repmat(rho_f,1,nx); tmp=tmp';
fid=fopen('prem_psv.den','wb');
fwrite(fid,tmp,'single');
fclose(fid);


fid=fopen('fd_par','w');
fprintf(fid,'toppad=%d\n',toppad);
fprintf(fid,'h=%f\n',h);
fprintf(fid,'nx=%d\n',nx);
fprintf(fid,'nz=%d\n',nz+toppad);
fprintf(fid,'zs=%d\n',ievdp+toppad);
fclose(fid);
%%
exit
