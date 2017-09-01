function dis2boundaryweight
    dirname='../pointnet.pytorch/real/data/real_half/2/';
    files=dir([dirname 'points/*.pts']);
    for i =100:size(files,1)
        v=dlmread([dirname 'points/' files(i).name]);
        f=dlmread([dirname 'faces/' strrep(files(i).name,'.pts','.fcs')]);
        f=f+1;
%        v(:,3)=-v(:,3);
        [w,vv]=addone(v,f,3,0.2,1);
%          dlmwrite([dirname 'weight/' strrep(files(i).name,'.pts','.wt')],w(:,1),' ');
%          dlmwrite([dirname 'alignedpoints/' (files(i).name)],vv,' ');
         display(i)
    end


%         f1=f1-1;e1=e1-1;
%         dlmwrite(['../pointnet.pytorch/real/data/real_half/1/points/1k_' filename(1:end-4) '__' num2str(cut_count) '1.pts'],v1,' ');
%         dlmwrite(['../pointnet.pytorch/real/data/real_half/1/faces/1k_' filename(1:end-4) '__' num2str(cut_count) '1.fcs'],f1,' ');
%         dlmwrite(['../pointnet.pytorch/real/data/real_half/1/edges/1k_' filename(1:end-4) '__' num2str(cut_count) '1.eds'],e1,' ');
%         f2=f2-1;e2=e2-1;
%         dlmwrite(['../pointnet.pytorch/real/data/real_half/2/points/1k_' filename(1:end-4) '__'  num2str(cut_count) '2.pts'],v2,' ');
%         dlmwrite(['../pointnet.pytorch/real/data/real_half/2/faces/1k_' filename(1:end-4) '__' num2str(cut_count) '2.fcs'],f2,' ');
%         dlmwrite(['../pointnet.pytorch/real/data/real_half/2/edges/1k_' filename(1:end-4) '__' num2str(cut_count) '2.eds'],e2,' ');


end

function [weight,v]=addone(v,f,xyz,initial,Align)
% halve a 3d mesh by x/y/z coordinate
%xyz= 1/2/3
%initial: a real number in the range of (0,1)
[e] = exterior_edges(f);
if (Align)
    v(e(:),3)=ones(length(e(:)),1)*mean(v(e(:),3));
end
[max_xyz,~]=max(v(:,xyz));
[vr,~]=find(v(:,xyz)>max_xyz-initial*(max(v(:,xyz))-min(v(:,xyz))));
weight=zeros(size(f));
want(1,1,:) =vr;
indexToDesiredRows = any(any(bsxfun(@eq,f,want),2),3);
rownum=find(indexToDesiredRows);
addedfaces=f(rownum,:);
weight(rownum,:)=weight(rownum,:)+1;
while size(addedfaces,1)<size(f,1)
% cutting line can be smoother
    Fr=zeros(0);
    for j=1:size(addedfaces,1)
            addface=addedfaces(j,:);
            for i=1:3
                [fr,~]=find(f==addface(i));
                Fr=[Fr;fr];
            end
    end
    Fr=unique(Fr,'rows');
    weight(Fr,:)=weight(Fr,:)+1;    
    newfaces=f(Fr,:);
%    end
    addedfaces=newfaces;
    addedfaces=unique(addedfaces,'rows');
end




%% plot
figure(2);
trisurf(f,v(:,1),v(:,2),v(:,3),'FaceVertexCData',weight/max(max(weight)));

end

