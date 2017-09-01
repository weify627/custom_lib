function plot_result_3d_training(dirname,time,pre)
%close all



%close all
    dirname='../pointnet.pytorch/half_real/results/20170802_172946/';
    pre='real_weighted0_summargin0barrier1area_cos10';
    dirname='../pointnet.pytorch/half_real/results/20170802_173311/';
    pre='real_weighted0_summargin0barrier1area_cos100000';
    %the three losses don't decrease but the picture looks better, what has
    %happened?
dirname='../pointnet.pytorch/half_real/results/20170802_231636/';
pre='real_20weightedv21000000000_summargin0barrier2000area_cos100000';

pre='res_norm_20weightedv21000000000_summargin0wbarrier2000area_cos100000';
dirname='../pointnet.pytorch/half_real/results/20170803_110853/';

dirname='../pointnet.pytorch/half_real/results/20170803_110406/';
pre='norm_2000weightedv21000000_summargin0wbarrier20000area_cos100000';

 dirname='../pointnet.pytorch/half_real/results/20170803_170951/';
 pre='xyz21000000000.0boundary20barrier20002000.0boundary0cos100000';

  dirname='../pointnet.pytorch/siamese/results/20170804_015620/';
 pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';
 
 dirname='../pointnet.pytorch/half_real/results/20170806_122401/';
 pre='signfirst_res0xyz41000000000boundary15_barrier2000002000boundary0cos100000';
%   dirname='../pointnet.pytorch/siamese/results/20170806_115127/';
%  pre='reimplementtest_res0xyz21000000000boundary20_barrier200000.02000boundary0cos100000sia100000';
 

 
%    dirname='../pointnet.pytorch/siamese/results/20170806_150717/';
%  pre='reimplementtest_res0xyz21000000000xyz20_barrier200000.02000boundary20cos100000sia100000';
%  
%  dirname='../pointnet.pytorch/siamese/results/20170806_152038/';
%    pre='reimplementtest_res0xyz210000000boundary20_barrier20000000.02000boundary20cos10000sia5000000';
  dirname='../pointnet.pytorch/half_real/results/20170806_155538/';
  pre='signfirst0xyz41000000000boundary15_barrier2000000boundary0cos100000';
  
    dirname='../pointnet.pytorch/half_real/results/20170806_164056/';
    pre='signfirst0xyz41000000000boundary15_barrier2000002000000000boundary0cos100000';
    
 %% change to margin loss (boundary weight  :bad
    dirname='../pointnet.pytorch/half_real/results/20170806_171106/';
    pre='signfirst0xyz41000000000boundary15_margin2000002000000000000000boundary0cos100000';
    
  %% margin loss (not boundary weight :bad
 dirname='../pointnet.pytorch/half_real/results/20170806_172224/';
pre='signfirst0xyz41000000000boundary15_margin2000002000000000000000N0cos100000'; 
  %% fine-tune with margin on previous best model
dirname='../pointnet.pytorch/half_real/results/20170806_220515/';
pre='ft_signfirst0xyz21000000000boundary20_margin2000002000000000boundary20cos100000';

dirname='../pointnet.pytorch/siamese/results/20170806_222607/';
pre='res0xyz61000000000xyz20_barrier200000.02000boundary20cos10000sia100000';
 

dirname='../pointnet.pytorch/half_real/results/20170806_231442/';
pre='ft_signfirst0xyz31000000000boundary20_margin2000002000000000000boundary20cos100000';
%    dirname='../pointnet.pytorch/half_real/results/20170805_230549/';
%  pre='signfirst_res0xyz41000000000boundary15_barrier2000002000boundary0cos100000';
%    dirname='../pointnet.pytorch/siamese/results/20170805_224603/';
%  pre='test_res0dis2boundary(0.6, 1.5)1000000000boundary15_barrier50000000.02000boundary0cos100000sia1000000';
%     dirname='../pointnet.pytorch/siamese/results/20170805_133900/';
%  pre='test_res0dis2boundary41000000000boundary15_barrier500000002000boundary0cos100000sia1000000';
% 
% dirname='../pointnet.pytorch/siamese/results/20170805_215519/';
% pre='test_res0dis2boundary(0.6, 1.5)1000000000boundary15_barrier50000000.02000boundary0cos100000sia1000000';
dirname='../pointnet.pytorch/siamese/results/20170807_174251/';
pre='res0xyz51000000000xyz20_barrier200000.02000boundary20cos100000sia100000';
% (dirname='../pointnet.pytorch/siamese/results/20170807_183706/';
% pre='res0xyz51000000000xyz20_margin200000.0200000000000boundary20cos100000sia100000';)
% 2*10^q,q=3 is too large; q=1 is good: whether or not normalize the
% s_boundary does not make much difference
% dirname='../pointnet.pytorch/siamese/results/20170807_190526/';
% pre='normpred_ptsres0xyz51000000000xyz20_margin200000.020boundary20cos100000sia100000';

%% harmonic
dirname='../pointnet.pytorch/harmonic/results/20170810_002756/';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';

%% harmonic +res area should be per-face scaled

% the followingfour trivial?
%overoverwelming margin+not weighted average for margin (left down)
dirname='../pointnet.pytorch/harmonic/results/20170812_003706/';
pre='stn1har100000_signfirst0xyz210000boundary20_margin200000.020000000000000N20cos100000';
%overwelming cos(right down)
dirname='../pointnet.pytorch/harmonic/results/20170811_235955/';
pre='stn1har_signfirst0xyz21000boundary20_margin200000.0200000boundary20cos100000';

%overwelming margin+har(right up)
dirname='../pointnet.pytorch/harmonic/results/20170812_001332/';
pre='stn1har_signfirst0xyz210000boundary20_margin200000.020000000000boundary20cos100000';
%overwelming margin(left up)
dirname='../pointnet.pytorch/harmonic/results/20170812_000854/';
pre='stn1har_signfirst0xyz210000boundary20_margin200000.020000000000boundary20cos100000';



%false:margin loss too large
dirname='../pointnet.pytorch/harmonic/results/20170812_135839/';
pre='res0stn1har100000_signfirst0xyz6100boundary20_margin2.02000000000000N20cos1';




%har+cos a bit small?
dirname='../pointnet.pytorch/harmonic/results/20170812_142831/';
pre='res0stn1har100000_signfirst0N00N20_margin200000.00N20cos100';

%har+area    +a bit margin or less area loss?
dirname='../pointnet.pytorch/harmonic/results/20170812_141649/';
pre='res0stn1har100000_signfirst0N010000N20_margin200000.00N20cos0';
%har+margin   
dirname='../pointnet.pytorch/harmonic/results/20170812_153434/';
pre='res0stn1har100000_signfirst0N00N20_margin200000.0200000N20cos0';
dirname='../pointnet.pytorch/harmonic/results/20170813_202654/';
pre='res0stn1har100000_signfirst0N0[1, 0]har0.001_margin200000.00N20cos0';

pre='res0stn1har100000_signfirst0N0[10000.0, 0]N0.001_margin200000.02000har0.01cos100';
bigtime='20170813_223735/';
pre='res0stn1har100000_signfirst0N0[10000.0, 0]N0.001_margin200000.00har0.01cos0';
bigtime='20170814_000548';
pre='res0stn1har100000_signfirst0N0[10000.0, 0]N0.001_margin200000.00har0.01cos0';
bigtime='20170814_003617';
%% half3400
pre='2mapres0stn1har100000000N0[10000000.0, 0.0]N0.001_margin200000.0200000har0.01cos100000b10000000';
bigtime='20170817_000702';

pre='2mapres0stn1har100000N0[10000.0, 0.0]N0.001_margin200000.0200har0.01cos100b10000';
bigtime='20170817_010540';
%commenting two modifications for /2/ in getitem is wrong
pre='2mapres0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.00har0.01cos0b0';
bigtime='20170817_003932';
% the following two are equivalant: 1st only comenting the z=-z,wrong
% flipped face
pre='2mapres0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.00har0.01cos0b0';
bigtime='20170817_003912';
%2nd add abs, right flipped face=>if margin loss is added, the 1st will be
%wrong
pre='2mapres0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.00har0.01cos0b0';
bigtime='20170817_012707';
% first fairly good combined
pre='2mapres0stn1har100000000N0[10000000.0, 0.0]N0.001_margin200000.0200000har0.01cos100000b10000000';
bigtime='20170816_132936';
pre='2mapl1res0stn1har0N0[10000000.0, 10000000.0]N0.001_margin200000.02e+13N20000cos100000b0';
bigtime='20170820_111339';

bigtime='20170821_091307';
pre='2mapl1res0stn1har0N0[10000000.0, 10000000.0]N0.001_margin200000.020000000000.0N20000cos100000b0';
%to check the following
bigtime='20170821_131951';
pre='smooth103d100000000000har0N0[100000000.0, 1000.0]N0.001_margin200000.00har0.01cos0b0';
bigtime='20170821_131910';
pre='smooth13d100000000000har0N0[100000000.0, 1000.0]N0.001_margin200000.00har0.01cos0b0';

% for the following, second is slightly better the the first => perhaps bn=0 for
% the last gconv layer is better
bigtime='20170824_174229/';
pre='2mapl1res0stn1har0N0[10000000000.0, 0.0]N0.001_margin200000.00.0N20000cos100000b0';
bigtime='20170824_164804/';
pre='2mapl1res0stn1har0N0[10000000000.0, 0.0]N0.001_margin200000.00.0N20000cos100000b0';

dirname=['../pointnet.pytorch/harmonic/results/' bigtime '/'];
% pre='2mapl1res0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.00har0.01cos0b0';
% bigtime='20170818_154314';
% dirname=['../pointnet.pytorch/smooth/results/' bigtime '/'];
 smooth=0;
plotr=3;
plotc=3;
winsize=[600 350 2000 2000];
winsize=[600 350 1500 1500];
num_points=1868;%1840;%704; %500 ;%1030;
num_face=3625;%3400;%1270; %1000; %2052;
labels = cellstr( num2str([1:num_points]') );
for i=1:5:num_points
    labels{i}=' ';
end
for k=16:1:62 %1000 %64
    fig = figure;
    u = fig.Units;
    fig.OuterPosition=(winsize);  


    time=1;
    fidin=fopen([dirname pre '_input_' num2str(time) '.txt']);
    in=textscan(fidin,'%f64 %f64 %f64');
    in=[in{1,1},in{1,2},in{1,3}];
    fidface=fopen([dirname pre '_face_' num2str(time) '.txt']);
    face=textscan(fidface,'%f %f %f');
    face=[face{1,1},face{1,2},face{1,3}];
    a_ori=in((k-1)*num_points+1:k*num_points,:);
    tri=face((k-1)*num_face+1:k*num_face,:);  
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

    %% calculate area
    plotn=1;
    in_area=triangle_area(a,tri);
       subplot(plotr,plotc,plotn); 
    hold on, trisurf(tri,a(:,1),a(:,2),a(:,3),'FaceColor','interp','FaceAlpha',0.6), hold off %,'FaceColor','interp'
%      text(a(1:max(max(tri)),1),a(1:max(max(tri)),2),a(1:max(max(tri)),3), labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
%                               'HorizontalAlignment','right');
    title(['ori '  num2str(k) ]);
    [E] = exterior_edges(tri);
    ori_a=a;
    if(k>32)
        ori_a=-a;
    end
    patch('Faces', tri, 'Vertices', a, 'FaceColor', [1,1,1], 'EdgeColor', [.7,.7,.7]);
    hold on;
    s = a(E(:,1),:); e = a(E(:,2),:);
    line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
    hold off;
    axis equal; axis tight; %axis off; 
    for time =1:5:3000%560:585%484:510 %460:485 % %250:275%149:173 %124:148 %100:123 %75:99 %50:74 %1:24 %49 %49
        plotn=plotn+1;
        fidout=fopen([dirname pre '_output2_' num2str(time) '.txt']);
        if smooth==0
            out=textscan(fidout,'%f%f', 'TreatAsEmpty','*','EmptyValue',10);
            out=[out{1,1},out{1,2}];
        else
            out=textscan(fidout,'%f%f%f', 'TreatAsEmpty','*','EmptyValue',10);
            out=[out{1,1},out{1,2},out{1,3}];
        end
        %close all
        b_ori=out((k-1)*num_points+1:k*num_points,:);
        b=b_ori(1:max(max(tri)),:);    
        out_area=triangle_area(b,tri);
        ratio=in_area./out_area;
        [Y,I]=sort(ratio);
        colormap=zeros(size(tri));
        for j=1:size(Y,1)
            colormap(I(j),:)=[j/size(Y,1),0.5*j/size(Y,1),1-j/size(Y,1)];
        end
     %%   figure plot
       subplot(plotr,plotc,plotn);  
        %subplot(1,4,3);
        axis equal;axis tight;
        grid on;
        a=b;
        %z=[1:max(max(tri))]'/max(max(tri));
        if smooth==0
            if k>32
                hold on, trisurf(tri,a(1:max(max(tri)),1),a(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.3,'EdgeColor','none','EdgeAlpha',0.2), hold off %,
            else
                hold on, trisurf(tri,a(1:max(max(tri)),1),a(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.3,'EdgeColor','none','EdgeAlpha',0.2), hold off %,
                
            end
        else
            hold on, trisurf(tri,a(1:max(max(tri)),1),a(1:max(max(tri)),2),a(:,3),'FaceColor','interp','FaceAlpha',0.6), hold off %,
        end
        title([ num2str(k) ' val' num2str(time)  ]);          
        
%         hold on, trisurf(tri,a(1:max(max(tri)),1),a(1:max(max(tri)),2),z,'FaceColor','interp','FaceAlpha',0.3), hold off %,
%     %         text(a(1:max(max(tri)),1), a(1:max(max(tri)),2),z, labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
%     %                                  'HorizontalAlignment','right');
%             title([ num2str(k) ' val' num2str(time)  ]);   
%         patch('Faces', tri, 'Vertices', a, 'FaceColor', [1,1,1], 'EdgeColor', [.7,.7,.7]);
%         hold on;
        s = a(E(:,1),:); e = a(E(:,2),:);
        if smooth==1
            line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'],'Color', 'r', 'LineWidth', 1);
        else
            line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], 'Color', 'r', 'LineWidth', 1);set(gca,'Color','k');
        end
        hold off;
        axis equal; axis tight; %axis off; 


        if plotn==plotr*plotc 
            plotn=0;
            suptitle(strrep(pre,'_',' '));
%            break
           w = waitforbuttonpress;
            if w == 0
                disp('Button click')
                %print('/home/weify/Pictures/real/trained_by_real/half/BarPlot','-dpng')
               % saveas(gcf,['/home/weify/Pictures/real/trained_by_real/half/' pre num2str(k) '.png']);
            else
                disp('Key press')
            end
               fig = figure;

            u = fig.Units;
            fig.OuterPosition=(winsize); 
       end
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
%   s23=[vin3-vin2 zeros(50,1)];
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
% training: trying to correct some overlapping
% dirname='../pointnet.pytorch/half_real/results/20170802_173311/';
% time='20170802_180148';
% pre='real_weighted0_summargin0barrier1area_cos100000';

% next step: barrier+cos?
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

