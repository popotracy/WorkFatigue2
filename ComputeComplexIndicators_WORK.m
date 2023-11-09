clear all; close all; clc;
addpath \\10.89.24.15\q\IRSST_DavidsTea\Raw_Data
addpath '\\10.89.24.15\q\IRSST_DavidsTea\Elvige'
addpath '\\10.89.24.15\j\IRSST_Fatigue\Pointage_repetitif_EMG\Scripts_IRSSTFatigue\Matlab_Fonctions\Entropy_measures'
addpath '\\10.89.24.15\j\IRSST_Fatigue\Pointage_repetitif_EMG\Scripts_IRSSTFatigue\Matlab_Fonctions\multiscaleSampleEntropy'
path_save = '\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators';
% chemin data markers: cd 'Q:\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered'
addpath(genpath('\\10.89.24.15\j\IRSST_Fatigue\Pointage_repetitif_EMG\Scripts_IRSSTFatigue\Matlab_Fonctions\Functions'))

Subjects = {'P01' 'P02' 'P03' 'P04' 'P05' 'P06' 'P07' 'P08' 'P09' 'P10'...
    'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20'...
    'P21' 'P22' 'P24' 'P25' 'P26' 'P27' 'P28' 'P29' 'P30' 'P31'};

Muscles = {'deltant';'deltmed';'deltpost';...
    'biceps';'triceps';'uptrap';'medtrap'; 'inftrap';'dent'};

Trials= {'Trial1';'Trial2';'Trial3';'Trial4'; ...
    'Trial5';'Trial6';'Trial7'; 'Trial8';'Trial9'};

Indicateurs = {'ActivationLevel','ApEntropy','CorrelationDimension',...
    'EMG_Activity','Fractal_Higuchi','FuzzyEntropy','HurstExponent','LypExpnt',...
    'Mobility','MSentropy','RR_DET','TFR_MedianFreq','TFR_SpectralEntropy'};

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

tic
for iSubject=1:length(Subjects)
    disp(Subjects{iSubject})
    Activity_WK = [] ;
    Mobility_WK = [] ;
    ApEntropy_WK = [] ; SamplEntropy_WK = []; MSentropy_WK = []; FuzzyEntropy_WK = [];
    FSS_Higuchi_WK = [];
    DOM_WK = [] ;
    CorDim_WK = [];
    REC_WK = [] ; DET_WK = [] ;
    LypExpnt_WK = [];
    
    load (['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\EMG_clean\WorkingTask\' Subjects{iSubject} '_Work']);
    load (['\\10.89.24.15\q\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered\Working_Task\' 'MatrixMarkersWork_' Subjects{iSubject}]);
    
    for iTrial = 1%1:length(Trials)
        disp(Trials{iTrial})
        Matrix = MatrixMarkersWork{1, iTrial} ; Matrix(1) = 1 ;
       
        Matrix = Matrix*1000/60 ;
        for iMuscle = 1%:length(Muscles)
            Data = EMG.data.(Trials{iTrial}).(Muscles{iMuscle}) ;
            for i=1:length(Matrix)-1
                Data_Markers = Data(round(Matrix(i)):(Matrix(i+1))) ;
                if sum(isnan(Data_Markers))==0 && sum(Data_Markers)~=0
                    Data_Markers_zscore  =  zscore(Data_Markers,1) ;
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
                else
                    %% Activity
                    Activity_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    %% Mobility
                    Mobility_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    %% Entropy Measures
                    ApEntropy_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    FuzzyEntropy_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    MSentropy_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    SamplEntropy_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    %% Fractals Self-Similarity: Higuchi
                    FSS_Higuchi_WK.(Muscles{iMuscle})(i,iTrial)  = nan ;
                    %% Multi-fractal detrended fluctuation: Hurst exponent & DOM
                    DOM_WK.(Muscles{iMuscle})(i,iTrial)= nan ;
                    %% Correlation Measures
                    CorDim_WK.(Muscles{iMuscle})(i,iTrial)= nan ;
                    REC_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    DET_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                    %% Chaos measure
                    LypExpnt_WK.(Muscles{iMuscle})(i,iTrial) = nan ;
                end
            end
        end
    end

    toc
    %% Data Saving
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\TimeFreq_Indicators\Activity\WK\' Subjects{iSubject} '_Activity_WK'],'Activity_WK')
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\TimeFreq_Indicators\Mobility\WK\' Subjects{iSubject} '_Mobility_WK'],'Mobility_WK')
%     
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\Entropy\SamplEntropy\WK\' Subjects{iSubject} '_SamplEntropy_WK'],'SamplEntropy_WK')
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\Entropy\MSEntropy\WK\' Subjects{iSubject} '_MSentropy_WK'],'MSentropy_WK')
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\Entropy\FuzzyEnt\WK\' Subjects{iSubject} '_FuzzyEntropy_WK'],'FuzzyEntropy_WK')
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\Entropy\ApEntropy\WK\' Subjects{iSubject} '_ApEntropy_WK'],'ApEntropy_WK')
%     
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\Higuchi\WK\' Subjects{iSubject} '_FSS_Higuchi_WK'],'FSS_Higuchi_WK')
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\DOMultifractality\WK\' Subjects{iSubject} '_DOM_WK'],'DOM_WK')
  
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\RQA\REC\WK\' Subjects{iSubject} '_REC_WK'],'REC_WK')
%     save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\RQA\DET\WK\' Subjects{iSubject} '_DET_WK'],'DET_WK')
    

     %save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\Complex_Indicators\CorDim\WK\' Subjects{iSubject} '_CorDim_WK'],'CorDim_WK')
     % 
     %   clear('ActLevel_RPT','Activity_RPT',...
     %        'MedianFreq_RPT','SpectralEntropy_RPT','Mobility_RPT',...
     %        'ApEntropy_RPT','FuzzyEntropy_RPT','SamplEntropy_RPT','MSentropy_RPT',...
     %        'FSS_Higuchi_RPT','DOM_RPT',...
     %        'CorDim_RPT','REC_RPT','DET_RPT','LypExpnt_RPT')
end