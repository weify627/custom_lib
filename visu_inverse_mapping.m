function visu_inverse_mapping(v,f,V)
% v: n*3 3D vertec matrix
% V: n * 2   2D map  e.g., harmonic
% xy,yq : x, y coordinate of griddata
[xq,yq] = meshgrid(-1:.2:1, -1:.2:1);
v1 = griddata(V(:,1),V(:,2),v(:,1),xq,yq);
v2 = griddata(V(:,1),V(:,2),v(:,2),xq,yq);
v3 = griddata(V(:,1),V(:,2),v(:,3),xq,yq);
dis=zeros(11,11);
for i=1:10
    for j=1:11
    dis(i,j)=(v1(i,j)-v1(i+1,j))^2+(v2(i,j)-v2(i+1,j))^2+(v3(i,j)-v3(i+1,j))^2;
    end
end
pink=[243, 222, 187]/255;
yellow=[253 ,201 ,207]/255;
figure()
trisurf(f,V(:,1),V(:,2),zeros(size(V,1),1),'FaceColor',[0.68,0.92,1],'FaceAlpha',1);hold on
scatter(xq(:),yq(:),10,'red','*');%ones(length(xq(:)),1)*yellow
figure()
hold on
trisurf(f,v(:,1),v(:,2),v(:,3),'FaceColor',[0.68,0.92,1],'FaceAlpha',1);
scatter3(v1(:),v2(:),v3(:),20,'red','*');
%ones(length(xq(:)),1)*yellow