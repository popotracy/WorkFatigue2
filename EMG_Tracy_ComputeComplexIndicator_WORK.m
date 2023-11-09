%% We want to know if the EMG data show fatigue before versus after the fatiguing repetitive pointing task.

clear, close all, clc




addpath(genpath('C:\Users\fabie\OneDrive\Bureau\Pont_serveur\btk'))
addpath(genpath('C:\Users\fh16095\Desktop\UdeM\WorkFatigue\Matlab_Fonctions\Functions'))
addpath(genpath('C:\Users\fh16095\Desktop\UdeM\WorkFatigue\Matlab_Fonctions'))
% addpath(genpath('\\10.89.24.15\j\IRSST_Fatigue\Pointage_repetitif_EMG\Scripts_IRSSTFatigue\Matlab_Fonctions\Functions'))


Participants = {'P1';'P2';'P4';'P5';'P6';'P7';'P8';'P9';'P10';'P11';'P12';'P13';'P14';'P15';'P16';'P17';...
    'P19';'P20';'P21';'P22';'P23';'P24'} ; %  ;'P3' (no post c3d file) ; 'P18' (no EMG) there are some bad participants, to be confirmed with Elvige

% the name of the muscle changes
Muscles = {'DeltA_IM_EMG5';'DeltM_IM_EMG6';'DeltP_IM_EMG7';'Bi_IM_EMG11';...
    'Tri_IM_EMG12';'Dent1_IM_EMG1';'TrapInf_IM_EMG10';'TrapMed_IM_EMG9';'TrapSup_IM_EMG8'} ; % up participant 17?
Muscles = {'Sensor_5_IM_EMG5';'Sensor_6_IM_EMG6';'Sensor_7_IM_EMG7';'Sensor_11_IM_EMG11';...
    'Sensor_12_IM_EMG12';'Sensor_1_IM_EMG1';'Sensor_10_IM_EMG10';'Sensor_9_IM_EMG9';'Sensor_8_IM_EMG8'} ; %  participant 18 to 24

Taskname='Tacheprefat'

for iP = 1:length(Participants)
    acq = btkReadAcquisition(['\\10.89.24.15\j\IRSST_Fatigue\Pointage_repetitif\Data\' Participants{iP} '\Trial\'Taskname'.c3d']) ; % Tachepostfat
    Data = btkGetAnalogs(acq) ;
    FreqSamp = btkGetAnalogFrequency(acq) ;

    EMG = [] ;
    for iM = 1:size(Muscles,1)
        EMG(:,iM) = Data.(Muscles{iM}) ;
    end
    figure(1) ; plot(EMG);

    %% clean data : check with Elvige that it is the same as for her analysis.
    [b,a] = butter(2,2*[10 400]/FreqSamp) ; % Parametre du filtre BP 10-400 Hz
    EMGBP = filtfilt(b,a,EMG) ;
    EMGBS = EMGBP;

    % Remove baseline
    EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)

    % Reducing sample frequency to 1000 Hz
    Ratio = FreqSamp / 1000 ;
    EMGDS = interp1(1:length(EMGBL),EMGBL,linspace(1,length(EMGBL),length(EMGBL)/Ratio));
    figure(2) ; plot(EMGDS);

    FreqSamp = 1000 ;     % Q: Should we downsample the frequency to 1000 from 2000?
    %EMG.Muscles = Muscles;
    %EMG.Borg = Data.Borg;
    save(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\WorkTask\EMG_WorkFatigue_Clean\' Participants{iP} '_' Taskname],'EMG')
end

%% Indicators parameters
r = 0.25 ; n = 2 ; % entropy measures
tau = 3 ; % Time Scale factor (source Wang et al 2021) window width was kept at 1s (2000 data points) with 50% overlap
m = 2 ; %m: match point(s) for MSE
Kmax = 8 ; %maximum number of sub-series composed from the original for Hihuchi FD
factor = 5 ; % factor: number of scale factor for MSE
dim = 2 ; % The embedding dimension (same as m for MSE)
scale = [2^2:1:16] ; % vector of scales for MFDFA
q = -10:0.1:10 ; pas = diff(q) ;% q-order that weights the local variations
m = 2 ; % polynomial order for the detrending
I = 1 ; % RQA
FreqSamp = 1000 ; 

tic

for iP=1:length(Participants)
    disp(Participants{iP})

    Activity_WK = [] ;
    Mobility_WK = [] ;
    ApEntropy_WK = [] ; SamplEntropy_WK = []; MSentropy_WK = []; FuzzyEntropy_WK = [];
    FSS_Higuchi_WK = [];
    DOM_WK = [];
    CorDim_WK = [];
    REC_WK = []; DET_WK = [];
    LypExpnt_WK = [];

    load (['C:\Users\fh16095\Desktop\UdeM\WorkFatigue\EMG_WorkFatigue_Clean\' Participants{iP} '_' Taskname '.mat']); % Lecture fichiers
    for iTrial=1 % pre=1 post=2
        for iM = 1:length(Muscles)
            Signal = EMG(:,iM);
            Signal = Signal(1:FreqSamp*100);
            Signal_zscore  =  zscore(Signal,1) ;

            %% Activity
            Activity_WK.(Muscles{iM})(iTrial) = var(Signal) ;

            %% Mobility
            dt= 1/FreqSamp ;
            first_derivative=diff(Signal)./dt;
            first_derivative(length(Signal))=first_derivative(length(first_derivative));
            var_derivative=var(first_derivative);
            var_EMG=var(Signal);
            Mobility_WK.(Muscles{iM})(iTrial) =sqrt(var_derivative./var_EMG);
            %% Entropy Measures: 4 indicators
            ApEntropy_WK.(Muscles{iM})(iTrial) = ApEn(Signal_zscore,dim,r) ;
            FuzzyEntropy_WK.(Muscles{iM})(iTrial) = FuzzyEn(Signal_zscore,dim,r,n) ;
            MSentropy_WK.(Muscles{iM})(iTrial) = multiscaleSampleEntropySignal_zscore,m,r,tau) ;
            SamplEntropy_WK.(Muscles{iM})(iTrial) = SampEn(Signal_zscore,dim,r) ;
            %% Fractals Self-Similarity: Higuchi
            FSS_Higuchi_WK.(Muscles{iM})(iTrial)  = Higuchi_FD(Signal,Kmax) ;
            %% Multi-fractal detrended fluctuation: Hurst exponent & DOM
            Hq = MFDFA1(Signal,scale,q,m,0);
            differentielle = (diff(Hq)./pas);
            interpdiff = interp1(1:length(differentielle),differentielle,linspace(1,length(differentielle),length(Hq)));
            alpha = Hq+q.*interpdiff;
            DOM_WK.(Muscles{iM})(iTrial)= max(alpha)-min(alpha);
            %% Correlation Measures: Correlation Dimension & RQA
            CorDim_WK.(Muscles{iM})(iTrial)= correlationDimension(Signal);

            RP= RPplot(Signal,15,4,2,1); % CD: Recurrence Quantification Analysis (RQA)
            REC_DET = Recu_RQA(RP,I);
            REC_WK.(Muscles{iM})(iTrial) = REC_DET(1) ;
            DET_WK.(Muscles{iM})(iTrial) = REC_DET(2) ;





        end

    end

    %% Compare pre versus post results: STATS: t test

    ttest(Mobility_pre,Mobility_post)
