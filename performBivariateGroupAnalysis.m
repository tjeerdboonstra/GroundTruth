function performBivariateGroupAnalysis

warning off

% Parameters
NSUR = 10000;
P = 0.05;
TYPE = 'Spearman';

TARGET = 5:7;
FEATURES = 8:75;

% Read data from table
data = readtable('groundtruth_meandata.csv');

targets = table2array(data(:,TARGET));
target_labels = data.Properties.VariableNames(TARGET);

features = table2array(data(:,FEATURES));
feature_labels = data.Properties.VariableNames(FEATURES);

pid = data.pid;
num_targets = length(TARGET);
num_features = length(FEATURES);
num_observations = length(pid);

% Compute correlations
[r,p] = corr(targets,features,'type',TYPE);


% Perform permutation test to estimate corrected p-value
fprintf('\nPermutation testing\n')

sur_r = zeros(num_targets,num_features,NSUR);
min_sur = zeros(NSUR,1);
max_sur = zeros(NSUR,1);
for n = 1:NSUR
    if ~rem(n,NSUR/10)
        fprintf('.')
    end
    i = randperm(num_observations);
    sur_features = features(i,:); % permute features
     
    sur_r(:,:,n) = corr(targets,sur_features,'type',TYPE);
    
    min_sur(n) = min(reshape(sur_r(:,:,n),1,[])); % select the lowest and high correlation for threshold
    max_sur(n) = max(reshape(sur_r(:,:,n),1,[]));    
end

% threshold is based on distribution of minimum and maximum correlations
min_sur = sort(min_sur);
max_sur = sort(max_sur);

lbound = min_sur(round(P*NSUR/2));
ubound = max_sur(round((1-P/2)*NSUR));


% Plot survey completion
figure(1)
set(gcf,'units','centimeters','position',[0 20 8.9 6])

bins = [0.5:1:18.5];
counts = histcounts(data.assessments,bins);

subplot('position',[0.12 0.14 0.78 0.78])
bar(bins(2:end)-0.5,cumsum(counts/sum(counts)*100,'reverse'))
set(gca,'xlim',[0.3 18.7], 'ylim',[0 100],'fontsize', 8,'xtick',[1:18],'ytick', [0:20:100],'ygrid','on');
xlabel('assessments','fontsize',9)
ylabel('percentage','fontsize',9)
title('completed assessments','fontsize',9)

yyaxis right
bar(bins(2:end)-0.5,cumsum(counts,'reverse'))
set(gca,'xlim',[0.3 18.7], 'ylim',[0 38],'fontsize',8,'xtick',[1:18],'ycolor',[0 0 0]);
ylabel('n','fontsize',9)
box off


% Plot results permutation test
figure(2)
set(gcf,'units','centimeters','position',[0 40 8.9 9])

subplot('position',[0.12 0.6 0.84 0.33])
[sorted_r,i] = sort(r(:));
sur_r = sort(reshape(sur_r,num_targets*num_features,NSUR),2);
plot(sur_r(i,:),'color',[.7 .7 .7])
hold on
line([1 numel(r)],[lbound,lbound],'color','k','linestyle','--')
line([1 numel(r)],[ubound,ubound],'color','k','linestyle','--')

plot(sorted_r,'k','linewidth',2)
i = find(sorted_r>ubound | sorted_r<lbound);
plot(i,sorted_r(i),'r.','markersize',20)
hold off

set(gca,'fontsize',8,'xlim',[0,numel(r)],'ylim',[-0.75, 0.75])
box off
title('Distribution of correlation coefficients','fontsize',9)
xlabel('comparisons','fontsize',9)
ylabel('rho','fontsize',9)

subplot('position',[0.12 0.1 0.84 0.33])
bins = linspace(-1,1,150);
min_h = hist(min_sur,bins);
max_h = hist(max_sur,bins);
bar(bins,min_h,'facecolor','b')
hold on
bar(bins,max_h,'facecolor','b')
line([lbound,lbound],[0 1000],'color','k','linestyle','--')
line([ubound,ubound],[0 1000],'color','k','linestyle','--')
hold off

set(gca,'xlim',[-0.75,0.75],'ylim',[0 1000],'box','off','fontsize',8)
title('Upper and lower threshold','fontsize',9)
xlabel('rho','fontsize',9)
ylabel('count','fontsize',9)


% Write significant results to command line
fprintf(['\n\nLower bound: ', num2str(lbound,'%0.2f'), ', upper bound: ', num2str(ubound,'%0.2f'),'\n',])

[i j] = find(r<lbound | r>ubound);
for n = 1:length(i)
    if r(i,j) < 0
        p_corr = mean(r(i(n),j(n))>=min_sur)*2;
    else
        p_corr = mean(r(i(n),j(n))<=max_sur)*2;
    end
    fprintf([target_labels{i(n)},' - ', feature_labels{j(n)},', r = ', num2str(r(i(n),j(n)),'%0.2f'), ', p_uncorr = ', num2str(p(i(n),j(n)),'%0.5f'), ', p_corr = ', num2str(p_corr,'%0.4f'),'\n'])
end


% Write all results to file
results = table;
counter = 1;
for t = 1:3
    [temp, i] = sort(r(t,:),'descend');
    
    for n = 1:length(i)
        results.target{counter} = target_labels{t};
        results.feature{counter} = feature_labels{i(n)};
        results.rho(counter) = r(t,i(n));
        results.CI_low(counter) = ci_r(t,i(n),1);
        results.CI_high(counter) = ci_r(t,i(n),2);
        results.p_uncorrected(counter) = p(t,i(n));
        if r(t,i(n))<0
            results.p_corrected(counter) = min(mean(r(t,i(n))>=min_sur)*2,1);
        else
            results.p_corrected(counter) = min(mean(r(t,i(n))<=max_sur)*2,1);
        end
        
        counter = counter+1;
    end
end
writetable(results,'groundtruth_correlations.csv')
            
            
    
    

