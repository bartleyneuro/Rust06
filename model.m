clear all
%define model parameters
V1.prefdirections = 30:30:360;  %list of each V1 neuron's preferred directions
STIM.possdirections = 30:30:360;  %list of possible stimulus directions

V1.tuningwidth = 180*pi/180; %tuning width of V1 neurons
V1.untunednormfactor = .1; %'untuned' normalization factor
V1.selfnormfactor = .1; %'self' normalization factor

MT.scalingnonlin = 1; %MT scaling nonlinearity
MT.exponentnonlin = 2; %MT exponent nonlinearity

%make up some weights of V1 influence on MT response:

%use this line for pattern cell
%V1.influenceweight = exp(-(linspace(-1,1,length(V1.prefdirections)).^2)/.4)-.3;

%use this line for a component cell
V1.influenceweight = exp(-(linspace(-1,1,length(V1.prefdirections)).^2)/.01); 
V1.influenceweight(3) = .2;

%Plot the linear weights of V1 neurons
figure(1);
clf
plot(V1.prefdirections,V1.influenceweight,'bo-','MarkerFaceColor','b');
xlabel('Preferred orientation (deg)');
ylabel('Weight');
xlabel('Direction (deg)');
set(gca,'XTick',0:45:360);
set(gca,'XLim',[0,360]);

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

