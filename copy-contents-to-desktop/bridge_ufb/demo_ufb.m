clear all;
close all;

load bridge_features_example.mat;

% Demonstrate calculation of user feedback from features
Nr = size(features,1);
fbn = [];
for i = 1:Nr
    if (all(~isnan(features(i,:))))
        fbn = [fbn; features_ufb(features(i,:),wvec,ilog,ftmin,ftmax,fbmin,fbmax)];
    end
end

figure;
subplot(2,1,1);
hist(fbn,50);
xlabel('Score');
ylabel('Count');
subplot(2,1,2);
plot(fbn,'k.');
xlabel('Test index');
ylabel('Score');
axis tight;
