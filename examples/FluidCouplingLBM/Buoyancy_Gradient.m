%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Luc Sibille  luc.sibille@3sr-grenoble.fr
%         Franck Lomine  franck.lomine@insa-rennes.fr
%
%   Script to interpret the DEM-LBM buoyancy simulation performed with the script buoyancy.py
%   This script was written for the OZ ALERT school in Grenoble 2011
%   Script updated in 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

radFactor = 0.6;
rhoSolid = 3000;
gravity = 0.1;

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
fscanf(fid,'%s',[1,1]);     rhoFluid=fscanf(fid,'%f',[1,1]);
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


C=dx/LBMdt;
Cs2=C^2/3;


path1 = './lbm-nodes/'
path2 = './dem-bodies/'

fiVx=dir([path1 'Vx_*']);
fiRho=dir([path1 'Rho_*']);
fiBodies=dir([path2 'spheres_*']);

nb_file =length(fiVx)




for i=1:nb_file
    
    
    vxfile=[path1 fiVx(i).name]
    ux = importdata(vxfile);
    ux=ux';
    
    Vdarcy(i) = mean( ux(1,:) );
    
    
    rhofile=[path1 fiRho(i).name];
    rho = importdata(rhofile);
    rho=rho';
    
    DP(i) = Cs2*(rho(1,90)-rho(end,90));
    
    
    fid2 = fopen([path2 fiBodies(i).name],'r');
    [donnees2] = fscanf(fid2,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[17,inf]);
    
    id=donnees2(1,:);
    x=donnees2(2,:);
    y=donnees2(3,:);
    r=donnees2(5,:);
    
    fclose(fid2);
    clear donnees2
    
    xmin(i)=min(x);
    xmax(i)=max(x);
    
    if i==1 %computation of Terzaghi critical hydraulic gradient
        
        Stotal = 0.005*(max(x+r)-min(x-r));
        Ssolid = sum( pi*(r/radFactor).^2 );
        
        Swater = Stotal-Ssolid;
        
        Wwater = rhoFluid*gravity*Swater;
        Wsolid = rhoSolid*gravity*Ssolid;
        %Wsolid = rhoSolid*gravity*sum( pi*(r/radFactor).^2 );
        
        Wtotal = Wwater+Wsolid;
        
        HgradCrit = Wtotal/Stotal/(rhoFluid*gravity);
    end
    
end

L = xmax(1)-xmin(1);
Hgrad = DP/(rhoFluid*gravity*L);

figure(1)
hold on
%plot(DP,Vdarcy)
plot(Hgrad,Vdarcy,'b')
plot([HgradCrit HgradCrit], [0 max(Vdarcy)],'r')
xlabel('Hydraulic gradient')
ylabel('Darcy velocity')
legend('Darcy velocity','Terzaghi critical gradient')

figure(2)
hold on
% plot(DP,xmin,'b')
% plot(DP,xmax,'r')
plot(Hgrad,xmin,'b')
plot(Hgrad,xmax,'g')
plot([HgradCrit HgradCrit], [0 max(xmax)],'r')
xlabel('Hydraulic gradient')
ylabel('Max and min particle positions')
legend('Top particle position','Bottom particle position','Terzaghi critical gradient')

