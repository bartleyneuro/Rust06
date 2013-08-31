%Now we will make the stimuli used to generate the matrix, Model, which is the 
%response to each pairwise combination of gratings that make all possible 
%plaids (as in figure 4).
count = 0;
STIM.gratingcomb = zeros(length(STIM.possdirections)^2,length(STIM.possdirections));
for i=1:length(STIM.possdirections)
    for j=1:length(STIM.possdirections);
        count = count+1;
        STIM.gratingcomb(count,i) = STIM.gratingcomb(count,i)+1/6;
        STIM.gratingcomb(count,j) = STIM.gratingcomb(count,j)+1/6;
    end
end
 
 
%Implement the MT model of Rust et al. (2006) to get the model responses
%------------------------------------------------------
clear Model
%1. Calulate the tuning curves for each neuron using the 'von Mises'
%function.  'TuningMatrix' will be a matrix with the rows are different neurons, and
%the columns are the response of each neuron to each stimulus orientation.
[thetam,pn] = meshgrid(STIM.possdirections,V1.prefdirections);
TuningMatrix = exp(V1.tuningwidth*cos( (thetam-pn)*pi/180));
 
%2. normalize so the sum of each row is 1.
TuningMatrixprime = TuningMatrix./ repmat(sum(TuningMatrix,2),1,length(STIM.possdirections));
 
%now, loop through each stimulus to get the model's response.
%Each row of STIM.gratingcomb is a different stimulus, and contains a list of contrasts 
%for each of the directions.
 
for i=1:size(STIM.gratingcomb,1)
 
    %3. Calculate the linear response to the stimulus
    LinearResponse = sum(TuningMatrixprime.*repmat(STIM.gratingcomb(i,:),length(V1.prefdirections),1),2);
 
    %4. Calculate the 'untuned' normalization
    UntunedNorm = LinearResponse.^2 ./ (sum(LinearResponse.^2) + V1.untunednormfactor^2);
    
    %5. Calcualte the 'self' normalization
    SelfNorm = UntunedNorm./(UntunedNorm+V1.selfnormfactor);
 
    %6. Calculate the linear sum across V1 neurons
    LinearSumV1s = V1.influenceweight*SelfNorm;
    
    %LinearSumV1s = V1.influenceweight*LinearResponse;   %no normalization
    %LinearSumV1s = V1.influenceweight*UntunedNorm;     %just untuned normalization
    
    %7. Calculate the nonlinear output
    Model(i) = MT.scalingnonlin*max(LinearSumV1s,0).^MT.exponentnonlin;
     %Model(i) = LinearSumV1s; %The input without MT nonlinearity
end
%-------------------------------------------------
 
%reshape it into a square
Model = reshape(Model,length(STIM.possdirections),length(STIM.possdirections));
 
%Display the matrix like they do in figure 4
figure(2);
clf
imagesc(STIM.possdirections,STIM.possdirections,Model);
colormap(gray)
hold on
contour(STIM.possdirections,STIM.possdirections,Model,8,'-','LineWidth',2,'Color',[1,.5,.5]);
set(gca,'YDir','normal');
axis equal
axis tight
xlabel('Direction 1 (deg)');
ylabel('Direction 2 (deg)');
set(gca,'XTick',0:90:360);
set(gca,'YTick',0:90:360);
 
figure(3)
clf
surf(STIM.possdirections,STIM.possdirections,Model);
set(gca,'XLim',[0,360]);
set(gca,'YLim',[0,360]);
set(gca,'XTick',0:90:360);
set(gca,'YTick',0:90:360);
 
