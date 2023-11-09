clear all; close all; clc;
addpath \\10.89.24.15\q\IRSST_DavidsTea\Raw_Data
addpath '\\10.89.24.15\q\IRSST_DavidsTea\Elvige'
addpath '\\10.89.24.15\j\IRSST_Fatigue\Pointage_repetitif_EMG\Scripts_IRSSTFatigue\Matlab_Fonctions'

% chemin data markers: cd 'Q:\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered'
addpath(genpath('V:\IRSST_Fatigue\Pointage_repetitif_EMG\Scripts_IRSSTFatigue\Matlab_Fonctions\Functions'))

Subjects = {'P01' 'P02' 'P03' 'P04' 'P05' 'P06' 'P07' 'P08' 'P09' 'P10'...
    'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20'...
    'P21' 'P22' 'P24' 'P25' 'P26' 'P27' 'P28' 'P29' 'P30' 'P31'};

Muscles = {'deltant';'deltmed';'deltpost';...
    'biceps';'triceps';'uptrap';'medtrap'; 'inftrap';'dent'};

Trials= {'Trial1';'Trial2';'Trial3';'Trial4'; ...
    'Trial5';'Trial6';'Trial7'; 'Trial8';'Trial9'};


for iSubject=1:length(Subjects)
    disp(Subjects{iSubject})
    MedFreq_SP=[]; MedFreq_RPT=[]; MedFreq_Work=[];
    EntSpect_SP=[];    EntSpect_RPT=[];  EntSpect_Work=[];
    load (['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\EMG_clean\SP\' Subjects{iSubject} '_SP']); % Lecture fichiers
    load (['\\10.89.24.15\q\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered\Static_Pose\' 'MatrixMarkersSP_' Subjects{iSubject}]);
    
%     load (['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\EMG_clean\RPT\' Subjects{iSubject} '_RPT']);
%     load (['\\10.89.24.15\q\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered\RPT\' 'MatrixMarkersRPT_' Subjects{iSubject}]);
%     
%     load (['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\EMG_clean\WorkingTask\' Subjects{iSubject} '_Work']);
%     load (['\\10.89.24.15\q\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered\Working_Task\' 'MatrixMarkersWork_' Subjects{iSubject}]);
    for iTrial = 1:length(Trials)
         disp(Trials{iTrial})
        for iMuscle = 1:length(Muscles)         
            Data_SP = EMG.data.(Trials{iTrial}).(Muscles{iMuscle}) ;
%             Data_RPT = EMG.data.(Trials{iTrial}).(Muscles{iMuscle}) ;
%             Data_WK =  EMG.data.(Trials{iTrial}).(Muscles{iMuscle});             
            %%  
        if sum(EMG.data.(Trials{iTrial}).(Muscles{iMuscle}))~=0 & nansum(EMG.data.(Trials{iTrial}).(Muscles{iMuscle}))~=0      
            %% Median frequency and Spectral entropy
            % TFR
%             if sum(Data_RPT(:,iMuscle))~=0 & nansum(Data_RPT(:,iMuscle))~=0  
                ond="cmor8-1";
                f=10:1:400 ;
                Fs = 1000; 
                TimeScale=Fs*centfrq(ond)./f;
                dt = 1/Fs;
                [coef,Freq]=cwt(Data_SP,TimeScale,ond,'ExtendSignal',1); %changer signal dentree en fonction de la tache
                TFR=abs(coef);
                MedianFreq= Compute_Median_Frequency(TFR,f) ;
                EntSpect = Compute_Spectral_Entropy(TFR,f) ;

                % Decoupage par cycle
                Matrix = MatrixMarkersSP{1, iTrial} ; Matrix(1) = 1 ;
%                 Matrix = MatrixMarkersRPT{1, iTrial} ; Matrix(1) = 1 ;
%                 Matrix = MatrixMarkersWork{1, iTrial} ; Matrix(1) = 1 ;

                Matrix = Matrix*1000/60 ;
                %% Static Pose
                for i=1:length(Matrix)
                    if isnan(Matrix(i))==0
                        if isnan(MedianFreq(round(Matrix(i)))) ==0
                            FreqMed(i,1)= MedianFreq(round(Matrix(i)));
                            SpectEnt(i,1)= EntSpect(round(Matrix(i)));
                        else
                            FreqMed(i,1) = nan ;
                            SpectEnt(i,1) = nan ;
                        end
                    else
                        FreqMed(i,1) = nan ;
                        SpectEnt(i,1) = nan ;
                    end
                end
              
            FreqMed = mean(FreqMed);
            SpectEnt = mean(SpectEnt);
            MedFreq_SP.(Trials{iTrial})(:,iMuscle)= FreqMed;
            EntSpect_SP.(Trials{iTrial})(:,iMuscle)= SpectEnt ;
            
%             MedFreq_RPT.(Trials{iTrial})(:,iMuscle)= FreqMed ;
%             EntSpect_RPT.(Trials{iTrial})(:,iMuscle)= SpectEnt ;
            
%             MedFreq_Work.(Trials{iTrial})(:,iMuscle) = FreqMed ;
%             EntSpect_Work.(Trials{iTrial})(:,iMuscle) = SpectEnt ;
        end
        end 
    end
save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\TimeFreq_Indicators\MedianFrequency\SP\' Subjects{iSubject} '_MedFreq_SP'],'MedFreq_SP')
save(['\\10.89.24.15\q\IRSST_DavidsTea\Elvige\TimeFreq_Indicators\SpectralEntropy\SP\' Subjects{iSubject} '_EntSpect_SP'],'EntSpect_SP')
end
