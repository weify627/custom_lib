% dirname='/home/weify/Desktop/codes/pointnet.pytorch/real/data/real_half/1/harmonic'
% filelist=dir(dirname);
% close all
% 
% 
% collection={filelist(3).name(4:8)};
% cnt=1;
% count=ones(1);
% for i = 4:length(filelist)
%     if strcmp(collection(cnt),filelist(i).name(4:8))
%         count(cnt)=count(cnt)+1;
%     else
%         collection=[collection,filelist(i).name(4:8)];
%         count=[count,1];
%         cnt=cnt+1;
%     end
%     
% end
% need=400./count;
% figure('visible','on');
% max_f=0;
% max_v=0;
% for i = 3:length(filelist)
%     v=dlmread([dirname '/' filelist(i).name]);
%     f=dlmread([strrep(dirname,'points','faces') '/' strrep(filelist(i).name,'.pts','.fcs')]);
%     f=f+1;   
%     a = doublearea(v,f);
%      if any(abs(a(:))<1e-8)==1
%         continue;
%      end
%     findex=find(strcmp(collection,filelist(i).name(4:8)));
%      for j=1:floor(need(findex))
%         v=rotate_3d(v); 
% 
% 
% %         figure(1);
% %         trisurf(f,v(:,1),v(:,2),v(:,3),'FaceColor','interp');
%         tri=f;
%             [E] = exterior_edges(tri);
%         s = v(E(:,1),:); e = v(E(:,2),:);
%         line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
%         hold off;
%         axis equal; axis tight; %axis off; 
%         V=harmonic_mex(v,f);
% %         figure(2);
% %         trisurf(f,V(:,1),V(:,2),v(:,3),'FaceColor','interp');    
%         max_v=max(max_v,size(v,1));
%         max_f=max(max_f,size(f,1));
%         E=E+1;
%         dlmwrite([strrep(dirname,'points','rotate_non0face_points') '/' strrep(filelist(i).name,'.pts',[num2str(j) '.pts'])],v,' ');
%         dlmwrite([strrep(dirname,'points','edges') '/' strrep(filelist(i).name,'.pts',[num2str(j) '.eds'])],E,' ');    
%         dlmwrite([strrep(dirname,'points','harmonic') '/' strrep(filelist(i).name,'.pts',[num2str(j) '.pts'])],V(:,1:2),' ');
%      end
%      display(i);
% end

%clear
%close all
%T = readtable('mesh_stats_copied.txt', 'Delimiter','\t');

%% rotate half
%/home/weify/Desktop/codes/pointnet.pytorch/real/matlab
 LIBIGL_DIR='/home/weify/Desktop/codes/libigl/';
cnt=1;
dirname='../pointnet.pytorch/real/data/half3000/1/points/';
filelist=dir(dirname);
collection={filelist(3).name(4:8)};
cnt=1;
count=ones(1);
nf_t=20;
for i = 3:length(filelist)
    if strcmp(collection(cnt),filelist(i).name(4:8))
        count(cnt)=count(cnt)+1;
    else
        collection=[collection,filelist(i).name(4:8)];
        count=[count,1];
        cnt=cnt+1;
    end
    
end
max_p=0;
each_n=50;
counter=1;
need=each_n./count;
dirname='../pointnet.pytorch/real/data/half3000/';
files=dir([dirname '1/points/*.pts']);
files2=dir([dirname '2/points/*.pts']);
newdirname='../pointnet.pytorch/real/data/half3000/';
dis=zeros(size(files,1),1);
figure('visible','off');
for i =1144:size(files,1)
    findex=find(strcmp(collection,files(i).name(4:8)));
%     if count(findex)>10 && count(findex)<100
%         continue;
%     end
%     if dis(findex)==1
%         continue
%     end
 %    dis(findex)=1;
    display(count(findex))
    display(findex);
    display(files(i).name);
%     if need(findex)<1
%         flag=floor(count(findex)/each_n);
%         counter=counter+1;
%         if mod(counter,flag)~=0
%             continue;
%         end
%     end
    pos=find(files(i).name=='o',1,'last');
    v=dlmread([dirname '1/points/' files(i).name]);
%    f=dlmread([dirname '1/faces/' strrep(files(i).name,'.pts','.fcs')]);
    f=dlmread([dirname '1/faces/' files(i).name(1:pos-2) '.fcs']);
    
    f=f+1;   
    a = doublearea(v,f);
     if any(abs(a(:))<1e-8)==1
        continue;
    end   
    v2=dlmread([dirname '2/points/' files2(i).name]);
 %   f2=dlmread([dirname '2/faces/' strrep(files2(i).name,'.pts','.fcs')]);
    f2=dlmread([dirname '2/faces/' files(i).name(1:pos-2) '.fcs']);
    
    f2=f2+1;   
    a2 = doublearea(v2,f2);    
    if any(abs(a2(:))<1e-8)==1
        continue;
    end
    Vu=v;
    Fu=f;
    Vu2=v2;
    Fu2=f2;
%     [Vu,Fu] = loop_mex(v,f,1);  
%     [Vu2,Fu2] = loop_mex(v2,f2,0); 
    [E] = exterior_edges(Fu);        
    if size(unique(E(:)),1)~=size(E,1) continue; end 
     [Vu,Fu] = decimate_mex(Vu,Fu,nf_t);    
%     [Vu2,Fu2] = decimate_mex(Vu2,Fu2,nf_t);    
%     V=bsxfun(@minus, Vu, mean(Vu));
%     V = bsxfun(@rdivide, V, std(V,0,1));
%     V2=bsxfun(@minus, Vu2, mean(Vu2));
%     V2 = bsxfun(@rdivide, V2, std(V2,0,1));    
%     
%     V=harmonic_mex(V,Fu);
%     V2=harmonic_mex(V2,Fu2); 

    for j=1:0%min(floor(need(findex))+1,20)
        
        theta = 2*pi*rand(1);     
        Vu=rotate_2d(Vu,theta);
        Vu2=rotate_2d(Vu2,theta);    
%         figure(5);
        trisurf(Fu2,Vu2(:,1),Vu2(:,2),Vu2(:,3),'FaceColor','interp');     axis equal    
%             figure(6);
        trisurf(Fu,Vu(:,1),Vu(:,2),Vu(:,3),'FaceColor','interp');        axis equal
%         figure(3);
        trisurf(Fu2,V2(:,1),V2(:,2),Vu2(:,3),'FaceColor','interp');
%         figure(4);
        trisurf(Fu,V(:,1),V(:,2),Vu(:,3),'FaceColor','interp');
        v=v2;V=V2;f=f2;
%         figure(3);
        trisurf(f,v(:,1),v(:,2),v(:,3),'FaceColor','interp');
%         figure(4);
        trisurf(f,V2(:,1),V2(:,2),v2(:,3),'FaceColor','interp');
        cnt=cnt+1;
        max_p=max(max(max_p,size(Vu,1)),size(Vu2,1));
    %          dlmwrite([dirname 'weight/' strrep(files(i).name,'.pts','.wt')],w(:,1),' ');
%         dlmwrite([newdirname '1/points/' strrep((files(i).name), '.pts', ['ro' num2str(j) '.pts'])],Vu,' ');
%         dlmwrite([newdirname '2/points/' strrep((files(i).name), '.pts', ['ro' num2str(j) '.pts'])],Vu2,' ');    
    end

%     dlmwrite([newdirname '1/harmonic/' (files(i).name)],V,' ');
%     dlmwrite([newdirname '2/harmonic/' (files(i).name)],V2,' ');    
    [E] = exterior_edges(Fu);        
    [E2] = exterior_edges(Fu2);
%     [Vu] = laplacian_smooth(Vu,Fu,'cotan',unique(E(:)),0.07,'implicit');
%     figure(7);
%     trisurf(Fu,Vu(:,1),Vu(:,2),Vu(:,3),'FaceColor','interp');    axis equal
    W=cotmatrix(Vu,Fu);
    M=massmatrix(Vu,Fu,'voronoi');
    Minv=spdiags(diag(M).^(-1),0,size(Vu,1),size(Vu,1));
    L=Minv*W
    a=adjacency_matrix(Fu)
    
%     dlmwrite([newdirname '1/l/' (files(i).name)],full(L),' ');    
% 
%     W=cotmatrix(Vu2,Fu2);
%     M=massmatrix(Vu2,Fu2,'voronoi');
%     Minv=spdiags(diag(M).^(-1),0,size(Vu2,1),size(Vu2,1));
%     L=Minv*W;
%     dlmwrite([newdirname '2/l/' (files(i).name)],full(L),' ');    
    
    
%     E=E-1;
%     E2=E2-1;
%     dlmwrite([newdirname '1/edges/' strrep((files(i).name), '.pts','.eds')],E,' ');
%     dlmwrite([newdirname '2/edges/' strrep((files(i).name), '.pts','.eds')],E2,' ');          
%     Fu=Fu-1;
%     Fu2=Fu2-1;
%     dlmwrite([newdirname '1/faces/' strrep((files(i).name), '.pts','.fcs')],Fu,' ');
%     dlmwrite([newdirname '2/faces/' strrep((files(i).name), '.pts','.fcs')],Fu2,' ');  
%     display(i);
%     display(findex);
end



%% move file
% dirname='/home/weify/Desktop/codes/pointnet.pytorch/real/data/half3000/1/';
% target='/1/';
% filelist=dir([dirname 'points']);
% for i=3:100:length(filelist) 
%     movefile([dirname '/points/' filelist(i).name], [strrep(dirname, target, ['/test' target]) '/points']);
%     pos=find(filelist(i).name=='o',1,'last');
% %    copyfile([dirname '/harmonic/' filelist(i).name(1:pos-2) '.pts'], [strrep(dirname, target, ['/test' target]) '/harmonic']);   
%     copyfile([dirname '/faces/' filelist(i).name(1:pos-2) '.fcs'],[strrep(dirname, target, ['/test' target]) '/faces']);
%     copyfile([dirname 'edges/' filelist(i).name(1:pos-2) '.eds'], [strrep(dirname, target, ['/test' target]) '/edges']);
% %    copyfile([dirname 'l/' filelist(i).name(1:pos-2) '.pts'], [strrep(dirname, target, ['/test' target]) '/l']);   
% 
%     display(i);
% end

%% check normrow(lap*p)
    u=laplacian_smooth(v,f,'uniform',[],0.01)
    Vu=v;
    Fu=f;
    W=cotmatrix(Vu,Fu);
    M=massmatrix(Vu,Fu,'voronoi');
    Minv=spdiags(diag(M).^(-1),0,size(Vu,1),size(Vu,1));
    L=Minv*W;
u=V;
    W=cotmatrix(u,Fu);
    M=massmatrix(u,Fu,'voronoi');
    Minv=spdiags(diag(M).^(-1),0,size(u,1),size(u,1));
    Lu=Minv*W;

    n=normrow(L*v);
    nu=normrow(Lu*u);
    n-nu <0;