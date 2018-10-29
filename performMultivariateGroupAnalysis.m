function performMultivariateGroupAnalysis

% Parameters
NSUR = 10000;
NUM_FEAT = 4;

TARGET = [6,5,7];
FEATURES = 8:75;


% Read data from table
data = readtable('groundtruth_meandata.csv');

targets = table2array(data(:,TARGET));
target_labels = data.Properties.VariableNames(TARGET);

features = table2array(data(:,FEATURES));
feature_labels = data.Properties.VariableNames(FEATURES);

pid = data.pid;
npid = length(pid);
num_targets = length(TARGET);
num_features = length(FEATURES);


figure
set(gcf,'units','centimeters','position',[5 20 18 13.5])

% Perform PLS
for t = 1:num_targets
    fprintf([target_labels{t},'\n'])
    x = zscore(features);
    y = zscore(targets(:,t));
    
    % compute PLS
    [XL,YL] = plsregress(x,y);
    
    % bootstrap
    for n = 1:NSUR
        i = ceil(rand(npid,1)*npid);
        [tempx,tempy] = plsregress(x(i,:),y(i,:));
        [XLsur(:,:,n),T,INX] = reorderModes(tempx,XL);
    end
    
    % and select robust features  
    zXL = mean(XLsur,3)./std(XLsur,0,3);
    [m i] = sort(abs(zXL(:,1)),'descend');
    i = i(1:NUM_FEAT);
    
    % run full model
    [XL_full,YL_full,XS_full,YS_full,beta_full,PCTVAR_full,MSE_full] = plsregress(x,y,length(i),'CV',5,'mcreps',NSUR);
    
    
    % the reduced model
    x = x(:,i);
    [XL_red,YL_red,XS_red,YS_red,beta_red,PCTVAR_red,MSE_red] = plsregress(x,y,length(i),'CV',5,'mcreps',NSUR);
    
    
    % and the reduced model with one component
    [XL_red1,YL_red1,XS_red1,YS_red1,beta_red1,PCTVAR_red1,MSE_red1] = plsregress(x,y,1);
    
    
    % Store regression model
    model(t).target = target_labels{t};
    model(t).beta = zeros(num_features,1);
    model(t).beta(i) = beta_red1(2:end);
    model(t).exp_var = PCTVAR_red1(2,1);
    
    
    % plot MSE of both models
    subplot('position',[0.09 1.03-0.32*t, 0.27 0.25])
    bar([0:length(i)],[MSE_full(2,:);MSE_red(2,:)]')
    hold on
    plot([-0.5,length(i)+0.5],[MSE_full(2,1),MSE_full(2,1)],'k--')
    hold off
    box off
    set(gca,'xlim',[-0.5,length(i)+0.5],'ylim',[0.5 2],'fontsize',8)
    hy = ylabel(target_labels{t},'fontsize',10,'fontweight','bold');
    set(hy, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
    xlabel('number of components')
    if t==1
        title('Mean square error','fontsize',9)
        lh = legend({'full model','reduced model'},'location','northwest');
        set(lh,'box','off','fontsize',8)
    end
    
    % plot feature weights of both models
    subplot('position',[0.435 1.03-0.32*t, 0.27 0.25])
    bar(beta_red1(2:end,1))
    set(gca,'xlim',[0.3,length(i)+0.7],'xticklabels',feature_labels(i),'fontsize',8)
    box off
    if t==1
        title('coefficients','fontsize',9)
    end
    
    % compare estimated and measured scores
    Yfit = [ones(npid,1),x]*beta_red1;
    subplot('position',[0.78 1.03-0.32*t, 0.25/18*13.5 0.25])
    plot(Yfit,y,'o')
    box off
    set(gca,'fontsize',8,'xlim',[-2 2], 'ylim',[-2.6 2.6])
    xlabel('predicted score')
    ylabel('observed score')
    
    % write results to command line
    fprintf(['Delta MSE full: ', num2str(MSE_full(2,2)-MSE_full(2,1),'%1.2f'),'\n']);
    fprintf(['Delta MSE reduced: ', num2str(MSE_red(2,2)-MSE_red(2,1),'%1.2f'),'\n']);
    
    [R,P,RLO,RUP]=corrcoef(Yfit,y);
    fprintf(['Correlation predicted-observed: r = ', num2str(R(1,2),'%1.2f'),',range ',num2str(RLO(1,2),'%1.2f'),' - ', num2str(RUP(1,2),'%1.2f'),'\n\n'])
end

save('model','model')



