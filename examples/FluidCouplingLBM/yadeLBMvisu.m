%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Luc Sibille  luc.sibille@3sr-grenoble.fr
%         Franck Lomine  franck.lomine@insa-rennes.fr
%
%   Script to visualize the fluid flow simulated with LBM written for a practical work in the OZ ALERT school in Grenoble 2011.
%   Script updated in 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The only thing you have to do is choose below the iteration number for
%%%% which you want to plot the fluid flow. This iteration number should
%%%% correspond to one the saving iteration appearing in the same of the
%%%% files saved in ./dem-bodies or ./lbm-nodes

IterNumber=600; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vxfile=sprintf('lbm-nodes/Vx_%.10d',IterNumber);
vyfile=sprintf('lbm-nodes/Vy_%.10d',IterNumber);
ux = importdata(vxfile);
uy = importdata(vyfile);
%ux=ux';
%uy=uy';
u = sqrt(ux.*ux+uy.*uy);


rhofile=sprintf('lbm-nodes/Rho_%.10d',IterNumber);
rho = importdata(rhofile);


nodebdfile=sprintf('lbm-nodes/NodeBoundary_%.10d',IterNumber);
nodeBD = importdata(nodebdfile);
nodeBD=nodeBD';


spherefile=sprintf('dem-bodies/spheres_%.10d',IterNumber);

tmpMat = importdata(spherefile);
Sid=tmpMat(:,1);
Sx=tmpMat(:,2);
Sy=tmpMat(:,3); 
Sz=tmpMat(:,4);
Sr=tmpMat(:,5);
clear tmpMat



%----- reading of LBM.log file
logfile='LBM.log';
fid = fopen(logfile);
tline = fgetl(fid);
fscanf(fid,'%s %s %s',[3,1]);
fscanf(fid,'%s',[1,1]);     Lx0=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     Ly0=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     depth=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     thickness=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     NbBodies=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     dP=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     Nu=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     Rho=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     dx=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     Nx=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     Ny=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     Nz=fscanf(fid,'%f',[1,1]);
tline = fgetl(fid);
tline = fgetl(fid);
fscanf(fid,'%s',[1,1]);     tau=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     omega=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     IterMax=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     SaveMode=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     IterSave=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     SaveGridRatio=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     IntervalLBM=fscanf(fid,'%f',[1,1]);
fscanf(fid,'%s',[1,1]);     LBMConvergenceCriterion=fscanf(fid,'%f',[1,1]);
tline = fgetl(fid);
tline = fgetl(fid);
fscanf(fid,'%s %s',[2,1]);     LBMdt=fscanf(fid,'%f',[1,1]);
fclose(fid);
%-------------------------------


fig1=figure(1)
set(gcf, 'color', 'white');
clf
nodeBD=flipdim(nodeBD',1);
imagesc(nodeBD);
title('Obstacle Boundaries');
axis equal;

axis off


%---------------------------------------
%   VELOCITIES (colormap)
%--------------------------------------
fig2=figure(2)
set(gcf, 'color', 'white');
%set(gca,'ZDir','reverse')
clf;
hold on

imagesc(u);
title('Velocity colormap');

umax_current=max(max(u));
if(umax_current==0.) umax_current=umax_current+1.e-9;end;
caxis([0. umax_current]);
colorbar('location','southoutside');

theta = -pi:2*pi/40:pi;
for aa=1:length(Sx)
    
    fill(ceil(Sx(aa)/dx)+Sr(aa)/dx*cos(theta),ceil(Sy(aa)/dx)+Sr(aa)/dx*sin(theta),[1 1 1]);
end

hold off;
axis equal;
axis off;



fig3=figure(3)
rhomax_current=max(max(rho));
rhomin_current=min(min(rho));
set(gcf, 'color', 'white');
%set(gca,'ZDir','reverse')
clf;
hold on

%rho=flipdim(rho,1)
imagesc(rho);
title('Rho');
if(rhomax_current==rhomin_current) rhomax_current=rhomax_current+1.e-9;end;
caxis([rhomin_current rhomax_current]);
colorbar;

theta = -pi:2*pi/40:pi;
for aa=1:length(Sx)
    
    %fill(ceil(Sx(aa)/dx)+Sr(aa)/dx*cos(theta),ceil(Sy(aa)/dx)+Sr(aa)/dx*sin(theta),[1 1 1]);
    fill((Sx(aa)/dx)+1+ceil(Sr(aa)/dx)*cos(theta),(Sy(aa)/dx)+1+ceil(Sr(aa)/dx)*sin(theta),[1 1 1]);
end
 
hold off;
axis equal;
axis off;
