%% Response to gratings vs. 120 deg plaids

%Single gratings
STIM.gratingcomb= zeros(length(STIM.possdirections),length(STIM.possdirections));
for i=1:length(STIM.possdirections)
    STIM.gratingcomb(i,i) = 1/6;
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
    
    %LinearSumV1s = V1.influenceweight*LinearResponse;  %no normalization
    %LinearSumV1s = V1.influenceweight*UntunedNorm;    %just untuned normalization
    
    %7. Calculate the nonlinear output
    Model(i) = MT.scalingnonlin*max(LinearSumV1s,0).^MT.exponentnonlin;
     %Model(i) = LinearSumV1s; %The input without MT nonlinearity
end
%-------------------------------------------------

figure(4)
clf
h1=plot(STIM.possdirections,Model,'r-o','MarkerFaceColor','r');
hold on


%Component prediction (linear sum of shifted responses to single gratings)

CompPred=Model([3:12 1:2]) + Model([11:12 1:10]);
h2= plot(STIM.possdirections,CompPred,'go-','MarkerFaceColor','g');


%response 120 deg plaids
STIM.gratingcomb= zeros(length(STIM.possdirections),length(STIM.possdirections));
for i=1:length(STIM.possdirections)
    STIM.gratingcomb(i,mod(i-3,length(STIM.possdirections))+1) = 1/6;
    STIM.gratingcomb(i,mod(i+1,length(STIM.possdirections))+1) = 1/6;
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
    
    %LinearSumV1s = V1.influenceweight*LinearResponse;  %no normalization
    %LinearSumV1s = V1.influenceweight*UntunedNorm;    %just untuned normalization
    
    %7. Calculate the nonlinear output
    Model(i) = MT.scalingnonlin*max(LinearSumV1s,0).^MT.exponentnonlin;
     %Model(i) = LinearSumV1s; %The input without MT nonlinearity
end
%-------------------------------------------------


h3= plot(STIM.possdirections,Model,'bo-','MarkerFaceColor','b');

%legend([h1,h2,h3],{'Single grating','Predicted component response','Model response'});
%xlabel('Direction (deg)');
%set(gca,'XTick',0:45:360);
%set(gca,'XLim',[0,360]);
%ylim = get(gca,'YLim');
%set(gca,'YLim',[0,ylim(2)]);

%ylabel('Response');
