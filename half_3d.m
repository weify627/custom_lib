function [v1,v2,f1,f2,e1,e2]=half_3d(v,f,xyz,initial,Align)
% halve a 3d mesh by x/y/z coordinate
%xyz= 1/2/3
%initial: a real number in the range of (0,1)
[max_xyz,~]=max(v(:,xyz));
[vr,~]=find(v(:,xyz)>max_xyz-initial*(max(v(:,xyz))-min(v(:,xyz))));

want(1,1,:) =vr;
indexToDesiredRows = any(any(bsxfun(@eq,f,want),2),3);
rownum=find(indexToDesiredRows);
addedfaces=f(rownum,:);
while size(addedfaces,1)<0.48*size(f,1)
% cutting line can be smoother
    
    [vr,~]=find(v(:,xyz)>=min(v(addedfaces(:),xyz)));
    want2=zeros(1,1,length(vr));
    want2(1,1,:) =vr;
    indexToDesiredRows = any(any(bsxfun(@eq,f,want2),2),3);
    t1=f(:,1);t2=f(:,2);t3=f(:,3);
    newfaces=[t1(indexToDesiredRows),t2(indexToDesiredRows),t3(indexToDesiredRows)];
%     rownum=find(indexToDesiredRows);
%     newfaces=f(rownum,:);
    newfaces=unique(newfaces,'rows');
    if((size(newfaces,1)-0.5*size(f,1))>0.5*size(f,1)-size(addedfaces,1))
        break;
    end
%    if(size(setdiff(newfaces,addedfaces,'rows'),1)==0)
    diffaces=setdiff(addedfaces,newfaces,'rows');
    for j=1:size(diffaces,1)
            addface=diffaces(j,:);
            for i=1:3
                [fr,~]=find(f==addface(i));
                newfaces=[newfaces;f(fr,:)];
            end
    end
%    end
    addedfaces=newfaces;
    addedfaces=unique(addedfaces,'rows');
end
f1=addedfaces;
f2 = setdiff(f,f1,'rows');
%[E] = exterior_edges(f1);
[y,I1]=sort(f1(:));
[C1,~,IC1] = unique(y);
f1(I1)=IC1;
v1=v(C1,:);

[y,I2]=sort(f2(:));
[C2,~,IC2] = unique(y);
f2(I2)=IC2;
v2=v(C2,:);

[e1] = exterior_edges(f1);
[e2] = exterior_edges(f2);
if (Align)
    v1(e1(:),3)=ones(length(e1(:)),1)*mean(v1(e1(:),3));
    v2(e2(:),3)=ones(length(e2(:)),1)*mean(v2(e2(:),3));
end


%% plot
% figure(1);
% trisurf(f1,v1(:,1),v1(:,2),v1(:,3));
% labels = cellstr( num2str([1:size(v1,1)]') );
% text(v1(:,1),v1(:,2),v1(:,3), labels, 'VerticalAlignment','bottom', ...
%                                   'HorizontalAlignment','right');
%  figure(2);
% trisurf(f2,v2(:,1),v2(:,2),v2(:,3));
% labels = cellstr( num2str([1:size(v2,1)]') );
% text(v2(:,1),v2(:,2),v2(:,3), labels, 'VerticalAlignment','bottom', ...
%                                   'HorizontalAlignment','right');
    % v=v1;
% [E] = exterior_edges(f1);
% V=v1;
% patch('Faces', f1, 'Vertices', V, 'FaceColor', [1,1,1], 'EdgeColor', [.7,.7,.7]);
% hold on;
% s = V(E(:,1),:); e = V(E(:,2),:);
% line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
% hold off;
% axis equal; axis tight;
% 
% 
% figure(2);
% [E] = exterior_edges(f2);
% V=v2;
% patch('Faces', f2, 'Vertices', V, 'FaceColor', [1,1,1], 'EdgeColor', [.7,.7,.7]);
% hold on;
% s = V(E(:,1),:); e = V(E(:,2),:);
% line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
% hold off;
% axis equal; axis tight;

end

