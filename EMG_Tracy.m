%% We want to know if the EMG data show fatigue before versus after the fatiguing repetitive pointing task.

clear, close all, clc

addpath(genpath('C:\Users\fabie\OneDrive\Bureau\Pont_serveur\btk'))


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

for iP=1:length(Participants)
    MedFreq_SP=[]; MedFreq_RPT=[]; MedFreq_Work=[];
    EntSpect_SP=[];    EntSpect_RPT=[];  EntSpect_Work=[];
    load (['C:\Users\fh16095\Desktop\UdeM\WorkFatigue\EMG_WorkFatigue_Clean\' Participants{iP} '_' Taskname '.mat']); % Lecture fichiers
    plot(EMG)
    for iM=1:length(Muscles)
        Data_SP = EMG(:,iM) ;
        if sum(Data_SP)~=0 & nansum(Data_SP)~=0
            ond="cmor8-1";
            f=10:1:400 ;
            Fs = 1000;
            TimeScale=Fs*centfrq(ond)./f;
            dt = 1/Fs;
            [coef,Freq]=cwt(Data_SP,TimeScale,ond,'ExtendSignal',1); %changer signal dentree en fonction de la tache
            TFR=abs(coef);
            MedianFreq= Compute_Median_Frequency(TFR,f) ;
            EntSpect = Compute_Spectral_Entropy(TFR,f) ;


            %% compute the EMG fatigue indicators for Tacheprefat and Tachepostfat
            % The scripts that compute the indicators of fatigue are here:

            % Mobility
            Mobility_pre(iP,iM) = diff(EMG)/(1/Freqsamp) ;
            % 14 other indicators
            %% Activity
            Activity_WK.(Muscles{iMuscle})(i,iTrial) = var(Data_Markers) ;
            %% Mobility
            dt= 1/EMG.FreqS ;
            first_derivative=diff(Data_Markers)./dt;
            first_derivative(length(Data_Markers))=first_derivative(length(first_derivative));
            var_derivative=var(first_derivative);
            var_EMG=var(Data_Markers);
            Mobility_WK.(Muscles{iMuscle})(i,iTrial) =sqrt(var_derivative./var_EMG);
            %% Entropy Measures: 4 indicators
            ApEntropy_WK.(Muscles{iMuscle})(i,iTrial) = ApEn(Data_Markers_zscore,dim,r) ;
            FuzzyEntropy_WK.(Muscles{iMuscle})(i,iTrial) = FuzzyEn(Data_Markers_zscore,dim,r,n) ;
            MSentropy_WK.(Muscles{iMuscle})(i,iTrial) = multiscaleSampleEntropy(Data_Markers_zscore,m,r,tau) ;
            SamplEntropy_WK.(Muscles{iMuscle})(i,iTrial) = SampEn(Data_Markers_zscore,dim,r) ;
            %% Fractals Self-Similarity: Higuchi
            FSS_Higuchi_WK.(Muscles{iMuscle})(i,iTrial)  = Higuchi_FD(Data_Markers,Kmax) ;
            %% Multi-fractal detrended fluctuation: Hurst exponent & DOM
            Hq = MFDFA1(Data_Markers,scale,q,m,0);
            differentielle = (diff(Hq)./pas);
            interpdiff = interp1(1:length(differentielle),differentielle,linspace(1,length(differentielle),length(Hq)));
            alpha = Hq+q.*interpdiff;
            DOM_WK.(Muscles{iMuscle})(i,iTrial)= max(alpha)-min(alpha);
            %% Correlation Measures: Correlation Dimension & RQA
            CorDim_WK.(Muscles{iMuscle})(i,iTrial)= correlationDimension(Data_Markers);

            RP= RPplot(Data_Markers,15,4,2,1); % CD: Recurrence Quantification Analysis (RQA)
            REC_DET = Recu_RQA(RP,I);
            REC_WK.(Muscles{iMuscle})(i,iTrial) = REC_DET(1) ;
            DET_WK.(Muscles{iMuscle})(i,iTrial) = REC_DET(2) ;
            %% Chaos measure
            LypExpnt_WK.(Muscles{iMuscle})(i,iTrial) = lyaprosen(Data_Markers,0.1);




        end

    end

    %% Compare pre versus post results: STATS: t test

    ttest(Mobility_pre,Mobility_post)
