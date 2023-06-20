% https://www.mathworks.com/help/simbio/ug/pk-pd-modeling-and-simulation-to-guide-dosing-strategy-for-antibiotics.html

% Load model
sbioloadproject('AntibacterialPKPD.sbproj', 'm1') ;

%% Select dose objects in the model
doseNames = {'250 mg bid', '250 mg tid', '500 mg bid', '500 mg tid'};
for iDoseGrp = 1:length(doseNames)
    doseRegimens(iDoseGrp) = sbioselect(m1, 'Name', doseNames{iDoseGrp}) ;
end

% Setup
nPatients    = 1000         ; % Number of patients per dosage group
nDoseGrps    = 4            ; % Number of tested dosage regimens

mu    = @(m,v) log(m^2/sqrt(v+m^2));
sigma = @(m,v) sqrt(log(v/m^2+1));
m     = @(typicalValue) typicalValue;
v     = @(typicalValue,CV) typicalValue^2*CV^2;

% Patient demographics
rng('default');
Wt        = normrnd(51.6, 11.8, nPatients , nDoseGrps )   ; % units: kg
Age       = normrnd(71.8, 11.9, nPatients , nDoseGrps )   ; % units: years
Scr_mu    =    mu(m(0.78), v(0.78,0.328));
Scr_sigma = sigma(m(0.78), v(0.78,0.328));
Scr       = lognrnd(Scr_mu, Scr_sigma, nPatients , nDoseGrps )   ; % units: ml/minute
% Gender ratio
id        = 1:nPatients*nDoseGrps                            ;
idFemale  = randsample(id , round(0.26*nDoseGrps*nPatients)) ; % 26% Female

CrCL            = (140 - Age).*Wt./(Scr*72)                 ; % units: ml/minute
CrCL(idFemale)  = CrCL(idFemale)*0.85                       ; % multiply by 0.85 for females

Central_mu    =    mu(m(7.64), v(7.64,0.20));
Central_sigma = sigma(m(7.64), v(7.64,0.20));
k12_mu        =    mu(m(1.59), v(1.59,0.20));
k12_sigma     = sigma(m(1.59), v(1.59,0.20));
k21_mu        =    mu(m(2.26), v(2.26, 0.2));
k21_sigma     = sigma(m(2.26), v(2.26, 0.2));

Central = lognrnd(Central_mu , Central_sigma, nPatients , nDoseGrps); % units: liter
k12     = lognrnd(k12_mu, k12_sigma, nPatients , nDoseGrps)         ; % units: 1/hour
k21     = lognrnd(k21_mu, k21_sigma, nPatients , nDoseGrps)         ; % units: 1/hour

CL      = 1.07*CrCL + 45.6 + normrnd(0,22,nPatients,nDoseGrps); % units: ml/minute

k1_mu      =    mu(m(5.59e-5), v(5.59e-5, 0.2));
k1_sigma   = sigma(m(5.59e-5), v(5.59e-5, 0.2));
k2_mu      =    mu(m(0.0297) , v(0.0297, 0.2));
k2_sigma   = sigma(m(0.0297) , v(0.0297, 0.2));
Kmax_mu    =    mu(m(3.50)   , v(3.50, 0.159));
Kmax_sigma = sigma(m(3.50)   , v(3.50, 0.159));

k1   = lognrnd(k1_mu, k1_sigma, nPatients , nDoseGrps)     ; % units: 1/hour
k2   = lognrnd(k2_mu, k2_sigma, nPatients , nDoseGrps)     ; % units: 1/hour
Kmax = lognrnd(Kmax_mu, Kmax_sigma, nPatients , nDoseGrps) ; % units: 1/hour

% Discrete distribution of MIC values based on 71 P. aeruginosa strains
micValue  = [0.0625, 0.125, 0.25, 0.5 , 1 , 2 , 4 , 8 , 16 , 32 ] ;
micFreq   = [   5  ,    8 ,  9  , 14  , 7 , 8 , 9 , 5 , 2  , 4  ] ;

% Sample MIC values from a discrete distribution using randsample
MIC = nan(nPatients, nDoseGrps) ; % preallocate
for iDoseGrp = 1:nDoseGrps
    MIC(:, iDoseGrp) = randsample(micValue , nPatients, true , micFreq);
end

KC50 = exp(-1.91 + 0.898*log(MIC) + 1.06*randn(nPatients , nDoseGrps)) ; % units: microgram/milliliter


params      = {'Central', 'k12', 'k21', 'CL', 'k1', 'k2', 'Kmax', 'KC50'};

observables = {'[Bacterial Growth Model].Growing',...
               '[Bacterial Growth Model].Resting'};

tempdose = sbiodose('dose');
tempdose.Target = 'Central.Drug';
tempdose.AmountUnits = 'milligram';
tempdose.TimeUnits = 'hour';
tempdose.DurationParameterName = 'TDose';

simfunc  = createSimFunction(m1,params,observables,tempdose,'UseParallel',true);

phi = cell(1,nDoseGrps);
for i = 1:nDoseGrps
    phi{i} = [Central(:,i),k12(:,i),k21(:,i), ...
              CL(:,i), k1(:,i), k2(:,i), ...
              Kmax(:,i), KC50(:,i)];
end

if isempty(gcp)
    parpool;
end

tObs      = 0:24:336       ; % hour
nTPoints  = length(tObs)   ; % Number of sampling points

% Preallocate
cfu          = nan(nTPoints,nPatients);
log10CFU     = cell(1,nDoseGrps) ;

for i = 1:nDoseGrps
    disp(['Simulating group ', num2str(i),'...'])

    % Get the dose table directly from an existing dose object for each
    % dosing regimen.
    doseTable = getTable(doseRegimens(i));

    % Simulate
    simdata = simfunc(phi{i},[],doseTable,tObs);

    % Sum of growing and resting bacterial counts for each patient
    for j = 1:nPatients
       cfu(:,j) = sum(simdata(j).Data,2);
    end
    % Store log-transformed counts for each dose group.
    log10CFU{i} = log10(cfu);
end

% Save results
log10CFU_250bid = log10CFU{1} ;
log10CFU_250tid = log10CFU{2} ;
log10CFU_500bid = log10CFU{3} ;
log10CFU_500tid = log10CFU{4} ;

delete(gcp('nocreate'));

hax1(1) = subplot(2,2,1)
plotCFUCount(tObs, log10CFU_250bid, 'a. Dose 250 bid' )
hax1(2) = subplot(2,2,2)
plotCFUCount(tObs, log10CFU_250tid, 'b. Dose 250 tid' )
hax1(3) = subplot(2,2,3)
plotCFUCount(tObs, log10CFU_500bid, 'c. Dose 500 bid' )
hax1(4) = subplot(2,2,4)
plotCFUCount(tObs, log10CFU_500tid, 'd. Dose 500 tid' )

% Link subplot axes
linkaxes(hax1)
