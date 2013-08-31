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
