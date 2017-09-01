function plot_result_3d_siamese(dirname,time,pre)
%close all
num_points=704; %700; %500 ;%1030;
num_face=1270; %1270; %1000; %2052;
k=1;
winsize=[600 50 2000 1300];
winsize=[600 50 1100 570];
%% how to solve the matching: denser point sampling, fix the mapped boundary shape(first normalize, then constrain to square)


%test last model with aligned boundary
dirname='../pointnet.pytorch/siamese/results/20170804_233941/';
time='20170804_233943';
pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';

%test #1 model with aligned boundary
dirname='../pointnet.pytorch/siamese/results/20170804_235616/';
time='20170804_235617';
pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';
%=> siamese is better
% good example:
dirname='../pointnet.pytorch/siamese/results/20170804_015620/';
time='20170804_053155';
pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';



 dirname='../pointnet.pytorch/siamese/results/';
% %1filename='20170805_215519/test_res0dis2boundary(0.6, 1.5)1000000000boundary15_barrier50000000.02000boundary0cos100000sia1000000_face_20170805_220948.txt';
% %2filename='20170805_224603/test_res0dis2boundary(0.6, 1.5)1000000000boundary15_barrier50000000.02000boundary0cos100000sia1000000_face_20170806_001648.txt';
% %3filename='20170806_002809/test_res0dis2boundary(0.6, 1.5)1000000000boundary20_barrier200000.02000boundary0cos100000sia1000000_face_20170806_002811.txt';
% 1>3:1 flattened and well preserve in center
% 2 similar to 3
% %4filename='20170805_224918/test_res0dis2boundary(0.6, 1.5)1000000000boundary20_barrier200000.02000boundary0cos100000sia1000000_face_20170806_002123.txt';
% %5filename='20170806_004610/test_res0dis2boundary(0.3, 3)100000000boundary20_barrier200000.0200000boundary0cos100000sia100000_face_20170806_020741.txt';
% 1>5=>larger barrierpara is better; neither 0.6&1.5 nor 0.3&3 is good, because the former area-changing range is too small to flatten the surface while the 
% latter aqueeze the center too much so that the original shape is severely damaged 
% %6filename='20170806_004400/test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000_face_20170806_020632.txt';
% 6 similar to 1
% (not aligned:)
% %filename='20170804_015620/test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000_face_20170804_053155.txt';
%good:
%problem: boundary are still not well aligned, ellipsoid is not mapped well
filename='20170806_115127/reimplementtest_res0xyz21000000000boundary20_barrier200000.02000boundary0cos100000sia100000_face_20170806_131000.txt';

%% to look at to try larger area_basepara
%#2 xyz as the weight of area loss (and with larger weight is better. (best till now)
filename='20170806_150717/reimplementtest_res0xyz21000000000xyz20_barrier200000.02000boundary20cos100000sia100000_face_20170806_163743.txt';
filename='20170807_000929/res0xyz51000000000xyz20_barrier200000.02000boundary20cos100000sia100000_face_20170807_012306.txt';
%failed=>siamese boundary loss cannot be too big, otherwise the boundary will become
%a line, and half of the surface is folded
%filename='20170806_152038/reimplementtest_res0xyz210000000boundary20_barrier20000000.02000boundary20cos10000sia5000000_face_20170806_165151.txt';
%failed (only a little modification on #2)
%next step: how to avoid boundary being a line?
%filename='20170806_222607/res0xyz61000000000xyz20_barrier200000.02000boundary20cos10000sia100000_face_20170806_234750.txt';
%% previous codes were wrong in pred_area, the following are all good results
% filename='20170807_174251/res0xyz51000000000xyz20_barrier200000.02000boundary20cos100000sia100000_face_20170807_183908.txt';
% filename='20170807_190526/normpred_ptsres0xyz51000000000xyz20_margin200000.020boundary20cos100000sia100000_face_20170807_202609.txt';
 filename='20170807_190153/res0xyz51000000000xyz20_margin200000.020boundary20cos100000sia100000_face_20170807_202301.txt';
pre=filename(17:end-25);


% dirname='../pointnet.pytorch/siamese/results/20170804_200724/';
% time='20170804_200730';
% pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';
% %old version
% fidin=fopen([dirname pre '_input_' time '.txt']);
% fidface=fopen([dirname pre '_face_' time '.txt']);
% fidout=fopen([dirname pre '_output_' time '.txt']);

fidin=fopen(strrep([dirname filename],'face','input'));
fidface=fopen([dirname filename]);
fidout=fopen(strrep([dirname filename],'face','output'));

in=textscan(fidin,'%f64 %f64 %f64');
in=[in{1,1},in{1,2},in{1,3}];

face=textscan(fidface,'%f %f %f');
face=[face{1,1},face{1,2},face{1,3}];

a=textscan(fidout,'%f%f', 'TreatAsEmpty','*','EmptyValue',10);
a=[a{1,1},a{1,2}];
[r,c]=find(a==10);
[rr,~]=find(a(:,2)~=10);
last=1;
out=a;
for i = 1:size(r,1)
    out(r(i)+1:r(i)+250,:)=a(r(i)+1:r(i)+250,:)/100;
end
out=out(rr,:);
fig = figure;u = fig.Units;fig.OuterPosition=(winsize); 
   
plotr=2;plotc=3;
pos=-1;
for kk =1:10:1000

    pos=pos+1;
    k=floor(kk/32)*64+mod(kk,32);
    a_ori=in((k-1)*num_points+1:k*num_points,:);
    b_ori=out((k-1)*num_points+1:k*num_points,:);
    tri=face((k-1)*num_face+1:k*num_face,:);      
    plotone(a_ori,b_ori,tri,plotr,plotc,0,num_points,kk);
    
    k=k+32;
    a_ori=in((k-1)*num_points+1:k*num_points,:);
    b_ori=out((k-1)*num_points+1:k*num_points,:);
    tri=face((k-1)*num_face+1:k*num_face,:);   
    plotone(a_ori,b_ori,tri,plotr,plotc,1,num_points,kk);    
    suptitle(strrep(pre,'_',' '));
    w = waitforbuttonpress;
    if w == 0
        disp('Button click')
    else
        disp('Key press')
    end
    fig = figure;u = fig.Units;fig.OuterPosition=(winsize); 

end
end

function area=triangle_area(pts,tri)

if size(pts,2)==2
    pts=[pts,zeros(size(pts,1),1)];
end
vin1=pts(tri(:,1),:);
vin2=pts(tri(:,2),:);
vin3=pts(tri(:,3),:);
s12=vin2-vin1 ;
s13=vin3-vin1 ; 
crossp=cross(s12', s13');
area=0.5*sum(crossp.^2,1);
    
area=area';
end
%%%
%% cubicinterp
% results/test_input_20170713_102559.txt
% dirname='../results/';
% time='20170713_102559';
% pre='test';
%trained by 3d_ellipsoid
% dirname='../../3d_ellipsoid/results/20170713_134849/';
% time='1';
% dirname='../results/20170713_153202/';
% time='20170713_155618';
% pre='test2';
% pre='cubicinterptest';
% % add mse of area +marginarea+cos
%     0.35z:first good result
% % dirname='../../cubicinterp/results/20170713_225654/';
% % time='20170713_232848';
% % pre='test_area';
%     randhalf dataset :good result  flatter the surface from inside the
%     sllipsoid(conclusion from >5 images),zero input has the same output
%     for each mesh
% dirname='../results/20170714_094322/';
% time='20170714_101639';
% pre='randhalf_area';
%     1z:mostly failed  perhaps need to panish small area,need to find the
%     boundary between failure and success in terms of lambda in lambda*z
%dirname='../results/20170714_095204/';
% time='20170714_102320';
% pre='cubicinterp1z_area';
% trained by v1map inputtrans_noglobal_3d_ellipsoid model bad
% dirname='../../3d_ellipsoid/results/';
% time='20170714_093259';
% pre='randhalfcos_marginarea';
% observation: squash along y (y+ or y-)

%partialy failed, only those parallel( or slightly tilt) to the x-y plane are the most
%successful ones =>perhaps the input transform is not good enough
% dirname='../results/20170714_153813/';  
% time='20170714_163057';
% pre='0.6z_area';

% dirname='../results/20170714_131129/';  %failure 25,29(the hole boundary is the most distorted part,but without inputtransform,this 29 looks just as normal as other notransform failure case,
%perhaps its failure is due to inputtransform's failure)
%perhaps they have learned the ability of keeping the consistensy for the
%same hole shape with different hole directions.(bu observing the mapped hole boundary's point order)
% time='20170714_131131';
% pre='over_0.8_area';

%trained with over_0.8,but tested with 0.6z noinputtransform, partially
%developed
% dirname='../results/20170714_214444/';  
% time='20170714_214446';
% pre='over_0.8_area';

%notrans    squash from the front view&from x- to x+;successful cases are
%mostly from this direction=>the input transform module learns to rotate
%the surface so that the hole is on the left
% dirname='../results/20170714_234706/';  
% time='20170714_234709';
% pre='over_0.8_area';
%   o2i1 mean fail25,41,45,52  hard47, has fixed case 29 in the last successful
%   one
% dirname='../results/20170717_105349/';  
% time='20170717_105350';
% pre='over_0.8_weightedarea_mean';
%o2i1 sum
% dirname='../results/20170717_105554/';
% time='20170717_105556';
% pre='over_0.8_weightedarea_mean';

% dirname='../../real10k/results/20170717_154803/';
% time='20170717_154804';
% pre='over_0.8_weightedarea_sum_real10k';

% dirname='../../real10k/results/20170717_193958/';
% time='20170717_200147';
% pre='bumpy_weightedarea_sum';
% 
% dirname='../../real10k/results/20170720_113734/';
% time='20170720_113736';
% pre='bumpy_weightedarea_sum';

%% trained by real 1k
%noweightedarea  most fail
% dirname='../pointnet.pytorch/real/results/20170731_235109/';
% time='20170801_023401';
% pre='real_weightedarea_sum_margin';
%most fail
% dirname='../pointnet.pytorch/real/results/20170801_001635/';
% time='20170801_072949';
% pre='real_weighted1area_sum_margin_cos0';

% same mesh, different cut: little difference
% 2 slightly different mesh: greatly different
% dirname='../pointnet.pytorch/real/results/20170801_002636/';
% time='20170801_030044';
% pre='real_weighted0area_sum_margin_cos1';


%% half real
% overlapped
% dirname='../pointnet.pytorch/half_real/results/20170802_172124/';
% time='20170802_174631';
% pre='real_weighted0_summargin1barrier0area_cos10';

% => barrier better than margin loss (at least the former is trying to flatten it), but some tends to be a line
% dirname='../pointnet.pytorch/half_real/results/20170802_172946/';
% time='20170802_175836';
% pre='real_weighted0_summargin0barrier1area_cos10';

%cos larger cos better than not, but the easiest failed(cylinder worse than
%yuanzhui),human is good
% dirname='../pointnet.pytorch/half_real/results/20170802_173311/';
% time='20170802_180148';
% pre='real_weighted0_summargin0barrier1area_cos100000';

% add barrier loss:become one line ==>cosine loss not functioning
% mostly not flipped
% dirname='../pointnet.pytorch/half_real/results/20170802_213816/';
% time='20170802_234513';
% pre='real_20weighted1000000000_summargin0barrier2000area_cos100000';
%barrier v2 is better,at least not like a line, but most human case
%failed(folded)
% dirname='../pointnet.pytorch/half_real/results/20170802_231636/';
% time='20170803_005227';
% pre='real_20weightedv21000000000_summargin0barrier2000area_cos100000';


% wrong normalization method
% (dirname='../pointnet.pytorch/half_real/results/20170803_023135/';
% time='20170803_045419';
% pre='norm_20weightedv21000000000_summargin0barrier2000area_cos100000';
% %results//_face_
% 
% dirname='../pointnet.pytorch/half_real/results/20170803_025551/';
% time='20170803_051445';
% pre='norm_20000weightedv21000000_summargin0wbarrier2000area_cos100000';)

% next step: barrier+cos?


%+normalize

%not good :has kept the abs(area), but many flipped
% time='20170803_131055';
% dirname='../pointnet.pytorch/half_real/results/20170803_110406/';
% pre='norm_2000weightedv21000000_summargin0wbarrier20000area_cos100000';

%good eg 121&131(hand) no flipping on the boundary
% time='20170803_132339';
% dirname='../pointnet.pytorch/half_real/results/20170803_110806/';
% pre='norm_20weightedv21000000000_summargin0wbarrier2000area_cos100000';



%+residual
%better than the next one in that less flip at guaidian, but worse in that
%the area in guaidian is quite small(but on the other hand, this particular 
%area property on guaidian preserves the property of curvature somehow)
% eg 41&61 to compare these two model
%convergence is very fast
% dirname='../pointnet.pytorch/half_real/results/20170803_110853/';
% time='20170803_132414';
% pre='res_norm_20weightedv21000000000_summargin0wbarrier2000area_cos100000';

%xyz area (this one is actually with res(but forgot to include 'res' in its name))   Good! (best till now)
%failure 1:cannot let the face too big at guaidian, even for boundary  eg.two, three
% when flipping happens, decrease the area loss; greater difference for
% xyz;
% failure 2 overlap of two non-nerighboring parts  eg.630
% time='20170803_183844';
% dirname='../pointnet.pytorch/half_real/results/20170803_170951/';
% pre='xyz21000000000.0boundary20barrier20002000.0boundary0cos100000';
%xyz(no res)
%almost as good as the last one=>residual makes no difference,res only set
%the mapping's x-y coordinate among all possible rotations
%by zooming in, we can observe that the overlapped area is quite small
% dirname='../pointnet.pytorch/half_real/results/20170803_211504/';
% time='20170803_224312';
% pre='res0xyz21000000000boundary20_barrier20002000boundary0cos100000';

%failed
% dirname='../pointnet.pytorch/siamese/results/';
% filename='20170805_133900/test_res0dis2boundary41000000000boundary15_barrier500000002000boundary0cos100000sia1000000_face_20170805_170759.txt';
%% 3d_ellipsoid
% folded
% time='20170708_154521'
% pre='mse_area'
% squeeze the ellipsoid to a flat surface
% marginarea_cos_input_20170708_204814.txt

% dirname='../results/'
% time='20170709_120330'
% pre='mse_cos'
% dirname='../results/'
% time='20170708_144725'
% pre='mse_cos'

%%%
function plotone(a_ori,b_ori,tri,plotr,plotc,pos,num_points,kk)
    labels = cellstr( num2str([1:num_points]') );
    for i=1:num_points
        if mod(i,4)~=0
            labels{i}=' ';
        end
    end
    %% process face
    [ni1,~,~]=find(tri(:,1));
    [ni2,~,~]=find(tri(:,2));
    ni=unique([ni1;ni2]);
    tri=tri(ni,:)+1;
    a=a_ori(1:max(max(tri)),:);
    b=b_ori(1:max(max(tri)),:);
    %% calculate area
    in_area=triangle_area(a,tri);
    out_area=triangle_area(b,tri);
    ratio=in_area./out_area;
    [Y,I]=sort(ratio);
    lib=jet(size(tri,1));
    color_map=zeros(size(tri));
    color_map(I,:)=lib;
%     for j=1:size(Y,1)
%         colormap(I(j),:)=[j/size(Y,1),0.5*j/size(Y,1),1-j/size(Y,1)];
%     end
 %%   figure plot
   subplot(plotr,plotc,pos*plotc+1);
    axis equal;axis tight;
    hold on, trisurf(tri,a(:,1),a(:,2),a(:,3),'FaceColor','interp','FaceAlpha',0.45), hold off %,'FaceColor','interp'
%     text(a(1:max(max(tri)),1),a(1:max(max(tri)),2),a(1:max(max(tri)),3), labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
%                           'HorizontalAlignment','right');
    title(['3D mesh #'  num2str(kk) ]);
    [E] = exterior_edges(tri);

%% plot
    s = a(E(:,1),:); e = a(E(:,2),:);
    line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
    hold off;
    axis equal; axis tight; %axis off; 
    ori_a=a;
   subplot(plotr,plotc,pos*plotc+2);  
    axis equal;axis tight;
    grid on;
    a=b;
    z=[1:max(max(tri))]';
    hold on, trisurf(tri,a(1:max(max(tri)),1),a(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.5,'EdgeColor','interp','EdgeAlpha',0.2), hold off %,
         text(a(1:max(max(tri)),1), a(1:max(max(tri)),2),z, labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
                                  'HorizontalAlignment','right');
        title(['2D map #'  num2str(kk)  ]);          set(gca,'Color','k');
    %patch('Faces', tri, 'Vertices', a, 'FaceColor', [1,1,1], 'EdgeColor', [.7,.7,.7]);
    %hold on;
    s = a(E(:,1),:); e = a(E(:,2),:);
%     line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], 'Color', 'r', 'LineWidth', 3);
%     hold off;
    axis equal; axis tight; %axis off; 
   subplot(plotr,plotc,pos*plotc+3);
    axis equal;axis tight;
    grid on;
    a=b;
    z=[1:max(max(tri))]';
    %z=ones(size(a,1),1);
    hold on, trisurf(tri,a(1:max(max(tri)),1),a(1:max(max(tri)),2),flipud(z),'FaceVertexCData',color_map,'EdgeColor','none'), hold off %,'FaceColor','interp'  ,'FaceAlpha',0.6
    title(['per-face area change #'  num2str(kk)  ]);  %,colormap jet ,colorbar
end