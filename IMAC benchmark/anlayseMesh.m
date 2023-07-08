close all
clc
clearvars

load('elements.mat')
elements(:,1) = [];
load('nodes.mat')

figure
PlotMesh(nodes,elements,0);
nNodes = size(nodes,1);

indx0 = nodes(:,1) == 0;

nodesx0 = nodes(indx0,:);

ynodesProfile = nodesx0(:,2);
znodesProfile = nodesx0(:,3);

indRestriciton = ynodesProfile>0;
ynodesProfile = ynodesProfile(indRestriciton);
znodesProfile = znodesProfile(indRestriciton);

figure
plot(ynodesProfile,znodesProfile,'-*');
grid on
axis('equal')

[ynodesProfile,I] = sort(ynodesProfile);
znodesProfile = znodesProfile(I);

a = [ynodesProfile,znodesProfile];

figure
plot(ynodesProfile,znodesProfile,'-x');
grid on
axis('equal')



%%
indz0 = nodes(:,3) == 0;
nodesz0 = nodes(indz0,:);
figure
plot(nodes(indz0,1),nodes(indz0,2),'*')
axis('equal')
