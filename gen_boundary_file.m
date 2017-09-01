dirname='/home/weify/Desktop/codes/pointnet.pytorch/real/data/real/points'
filelist=dir(dirname);
close all


collection={filelist(3).name(4:8)};
cnt=1;
count=ones(1);
for i = 4:length(filelist)
    if strcmp(collection(cnt),filelist(i).name(4:8))
        count(cnt)=count(cnt)+1;
    else
        collection=[collection,filelist(i).name(4:8)];
        count=[count,1];
        cnt=cnt+1;
    end
    
end
need=300./count;
figure('visible','off');
max_f=0;
max_v=0;
for i = 8144:length(filelist)
    v=dlmread([dirname '/' filelist(i).name]);
    f=dlmread([strrep(dirname,'points','faces') '/' strrep(filelist(i).name,'.pts','.fcs')]);
    f=f+1;   
    a = doublearea(v,f);
     if any(abs(a(:))<1e-8)==1
        continue;
     end
    findex=find(strcmp(collection,filelist(i).name(4:8)));
     for j=1:floor(need(findex))
        v=rotate_3d(v); 


%         figure(1);
%         trisurf(f,v(:,1),v(:,2),v(:,3),'FaceColor','interp');
        tri=f;
            [E] = exterior_edges(tri);
        s = v(E(:,1),:); e = v(E(:,2),:);
        line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
        hold off;
        axis equal; axis tight; %axis off; 
        V=harmonic_mex(v,f);
%         figure(2);
%         trisurf(f,V(:,1),V(:,2),v(:,3),'FaceColor','interp');    
        max_v=max(max_v,size(v,1));
        max_f=max(max_f,size(f,1));
        E=E+1;
        dlmwrite([strrep(dirname,'points','rotate_non0face_points') '/' strrep(filelist(i).name,'.pts',[num2str(j) '.pts'])],v,' ');
        dlmwrite([strrep(dirname,'points','edges') '/' strrep(filelist(i).name,'.pts',[num2str(j) '.eds'])],E,' ');    
        dlmwrite([strrep(dirname,'points','harmonic') '/' strrep(filelist(i).name,'.pts',[num2str(j) '.pts'])],V(:,1:2),' ');
     end
     display(i);
end