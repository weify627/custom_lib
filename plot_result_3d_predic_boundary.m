function plot_result_3d_predic_boundary(dirname,time,pre)
%close all
num_points=550 %704; %700; %500 ;%1030;
num_face=1000 %1270; %1270; %1000; %2052;

labels = cellstr( num2str([1:num_points]') );
for i=1:num_points
    if mod(i,4)~=0
        labels{i}=' ';
    end
end
k=1;
winsize=[600 50 2000 1300];
winsize=[600 50 900 1300];
% 
% 
dirname='../pointnet.pytorch/half_real/results/20170806_231442/';
time='20170806_235646';
pre='ft_signfirst0xyz31000000000boundary20_margin2000002000000000000boundary20cos100000';

% comparing hte following two, finetuned is slightly better (flipped area is less)
% #1:time='20170803_183844';
% dirname='../pointnet.pytorch/half_real/results/20170803_170951/';
% pre='xyz21000000000.0boundary20barrier20002000.0boundary0cos100000';
% ft of #1
% dirname='../pointnet.pytorch/half_real/results/20170806_220515/';
% pre='ft_signfirst0xyz21000000000boundary20_margin2000002000000000boundary20cos100000';
% time='20170806_224344';

% area of faces in the center become larger, but why?
% dirname='../pointnet.pytorch/half_real/results/20170806_231442/';
% time='20170806_235646';
% pre='ft_signfirst0xyz31000000000boundary20_margin2000002000000000000boundary20cos100000';

%% to check: l_area.basepara=2
% dirname='../pointnet.pytorch/half_real/results/20170807_002104/';
% time='20170807_011729';
% pre='ft_signfirst0xyz21000000000boundary20_margin2000002000000000000boundary20cos100000';

%% to check: l_area.basepara=0.2(right bottom)
% dirname='../pointnet.pytorch/half_real/results/20170807_002128/';
% time='20170807_011821';
% pre='ft_signfirst0xyz01000000000boundary20_margin2000002000000000000boundary20cos100000';
% for the above two , don't know why all scaled down? ->10^(-3)
dirname='../pointnet.pytorch/harmonic/results/20170807_233409/';
time='20170808_000749';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
%% harmonic

dirname='../pointnet.pytorch/harmonic/results/20170810_102632/';
time='20170810_110114';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
% +26f (three and human are slightly better, but others are not)
dirname='../pointnet.pytorch/harmonic/results/20170810_234902/';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
time='20170811_024911';

dirname='../pointnet.pytorch/harmonic/results/20170811_015144/';
pre='stn1har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
time='20170811_023225';

%% cut mesh
dirname='../pointnet.pytorch/harmonic/results/20170812_190539/';
pre='res0stn1har100000_signfirst0N00N20_margin200000.02000N20cos0';
time='20170812_192438';
dirname='../pointnet.pytorch/harmonic/results/20170812_212041/';
pre='res0stn1har100000_signfirst0N0[10000.0, 1000.0]N20_margin200000.00N20cos0';
time='20170812_213353';  %% generally miss in pair, predict ~200,hit 2/3

fidhar=fopen([dirname pre '_har_' time '.txt']);
har=textscan(fidhar,'%f%f', 'TreatAsEmpty','*','EmptyValue',10);
har=[har{1,1},har{1,2}];
fidin=fopen([dirname pre '_input_' time '.txt']);
in=textscan(fidin,'%f64 %f64 %f64');
in=[in{1,1},in{1,2},in{1,3}];
fidface=fopen([dirname pre '_face_' time '.txt']);
face=textscan(fidface,'%f %f %f');
face=[face{1,1},face{1,2},face{1,3}];
fidout=fopen([dirname pre '_output_' time '.txt']);
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
   
plotr=5;plotc=2;
pos=-1;
for k =11:10:1000

    pos=pos+1;
    a_ori=in((k-1)*num_points+1:k*num_points,:);
    b_ori=out((k-1)*num_points+1:k*num_points,:);
    tri=face((k-1)*num_face+1:k*num_face,:);  
    b_har=har((k-1)*num_points+1:k*num_points,:);
    %% process face
    l=0;
    for j=100:num_face
        if tri(j,:)==zeros(1,3)
           l=j-1;
           break;
        end
    end
    if l==0
        l=num_face;
    end
    tri=tri(1:l,:);
    tri=tri+1;
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
   z=[1:max(max(tri))]';
    hold on, trisurf(tri,a(:,1),a(:,2),a(:,3),'FaceColor','interp','FaceAlpha',0.45), hold off %,'FaceColor','interp'
%     text(a(1:max(max(tri)),1),a(1:max(max(tri)),2),a(1:max(max(tri)),3), labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
%                           'HorizontalAlignment','right');
    title(['3D mesh '  num2str(k) ]);
    [E] = exterior_edges(tri);
    s = a(E(:,1),:); e = a(E(:,2),:);
    qq=length(s)
    line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
    hold off;
    axis equal; axis tight; %axis off; 
    ori_a=a;
%% plot
 subplot(plotr,plotc,pos*plotc+2);  
    axis equal;axis tight;
    grid on;
    [Y,I]=max(b(),[],2);
    vindex=find(I==2);
  %  E_predic1=a(vindex,:);
    hold on, trisurf(tri,a(:,1),a(:,2),a(:,3),'FaceColor','white','FaceAlpha',1)%, hold off 
                
    E_predic1=a(vindex(find(ismember(vindex,E))),:);
    missed=setdiff(E,vindex);
    false=setdiff(vindex,E);
    E_predic=a(missed,:);    
    singlemiss=size(E_predic,1)-size(unique(E_predic,'rows'),1);
    scatter3(E_predic(:,1),E_predic(:,2),E_predic(:,3),7,'red')
    E_predict=a(false,:);    
    scatter3(E_predict(:,1),E_predict(:,2),E_predict(:,3),7,'blue','*')
    scatter3(E_predic1(:,1),E_predic1(:,2),E_predic1(:,3),10,'yellow')
    title(['total:' num2str(length(s)) ', miss:'  num2str(length(missed)) '(' num2str(singlemiss) ' pairs), false:' num2str(length(false)) ]);
    if pos==plotr-1   %  fail25,41,45,52  hard47
        pos=-1;
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

%#1 xyz area (this one is actually with res(but forgot to include 'res' in its name))   Good! (best till now)
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


%assume the input faces are mostly the same size, differences among faces
%in the mapping will also indicate lots of info w.r.t. the area change of
%each face
% area of the highest faces are larger than #1
% dirname='../pointnet.pytorch/half_real/results/20170804_025604/';
% time='20170804_051902';
% pre='res0xyz41000000000boundary15_barrier20002000boundary0cos100000';

%compared to the last one, this one has larger infinite for BarrierArea (previous's basepara is 100000(not 2000))
%=> larger barrier loss (vs area loss is better)
%it preserves more original info of the points farthest to origin  e.g.
%four && larger (high face area)/(low face area) ratio
%less twirled e.g. case 111 (compared to last version)
%
%SSolution to do: boundary to curved ranking as input or supervision (perpoint)
% set the ranking fixed in a range of points
%SSolution: larger barrier loss is better
%dirname='../pointnet.pytorch/half_real/results/20170804_030107/';
% time='20170804_052355';
% pre='res0xyz41000000000boundary15_barrier200000002000boundary0cos100000';
% 
%till now typical failure: 201,one,two,three,four
% 
% dirname='../pointnet.pytorch/siamese/results/20170804_015620/';
% time='20170804_053155';
% pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';

%%cos1 vs cos100000 makes little difference
% dirname='../pointnet.pytorch/half_real/results/20170805_134937/';
% time='20170805_155020';
% pre='res0xyz41000000000boundary15_barrier200000002000boundary0cos1';

%the following two failed
% (dirname='../pointnet.pytorch/half_real/results/20170806_172224/';
% time='20170806_213507';
% pre='signfirst0xyz41000000000boundary15_margin2000002000000000000000N0cos100000';
% 
% dirname='../pointnet.pytorch/half_real/results/20170806_171106/';
% pre='signfirst0xyz41000000000boundary15_margin2000002000000000000000boundary0cos100000';
% time='20170806_213749';)

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

