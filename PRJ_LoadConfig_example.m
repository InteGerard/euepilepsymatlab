function [channels, se_alphas, se_bands, plv_alphas, plv_bands, downsample_rate, records_to_process] = EU_LoadConfig(patient_id)


se_alphas = [6 8 10 12 14];

se_bands = [ 0.5 4;
            4 8;
            8 13;
            13 30;
            30 60  ];
        
plv_alphas = [6 8 10 12 14];     

plv_bands = [ 0.5 4;
            4 8;
            8 12;
            12 25;
            25 45  ];
            
            
%% Configuration for each patient
if strcmp(patient_id, 'pat_FR_253')

    channels = { ...
        'HRB2' ; ...
        'HRC2' ; ...
        'HRB3' ; ...
        'HRA1' ; ...
        'HRA2' ; ...
        'HRA3' ; ...
        'HRB1' ; ...
        'HRC1' ; ...
        'HRC3' ; ...
        'HLC1' ; ...
        'HLB2' ; ...
        'HRA5' ; ...
        'HLA2' ; ...
        'HLA1' ; ...
        'HLA3' ; ...
        'HLA4'   ...
    };

     downsample_rate = 2;
     
     records_to_process = 1:999;
end   


%% Configuration for pat_FR_1096
if  strcmp(patient_id, 'pat_FR_1096')

    channels = { ...
        'HL2' ; ... %onsets
        'HL3' ; ...
        'HL4' ; ...
        'HL5' ; ...
        'HL6' ; ...
        'HL7' ; ...
        'HL8' ; ...
        'HL9' ; ... 
        'GA1' ; ... %early onsets
        'GB4' ; ...
        'GB5' ; ...
        'GB6' ; ...
        'GB7' ; ...
        'GB8' ; ...
        'GC1' ; ...
        'GC2' ; ...      
    };

     downsample_rate = 4;
     
     records_to_process = 1:9999;
end  

%% Configuration for pat_FR_970
if  strcmp(patient_id, 'pat_FR_970')

    channels = { ...
        'GF4'; ...
        'GG1'; ...
        'GG2'; ...
        'GG3'; ...
        'GG4'; ...
        'GH1'; ...
        'GH2'; ...
        'GH3'; ...
        'TBA1'; ...
        'TBA2'; ...
        'TBA4'; ...
        'TBB2'; ...
        'TL1'; ...
        'TL2'; ...
        'TL3'; ...
        'TL4'; ... 
    };

     downsample_rate = 1;
     
     %records_to_process = 25:182;
     records_to_process = 89:182;
end  
   

%% Configuration for pat_FR_548
if  strcmp(patient_id, 'pat_FR_548')

    channels = { ...
        'HL10'; ...
        'HL7'; ...
        'HL8'; ...
        'HL9'; ...
        'TBLA1'; ...
        'TBLA2'; ...
        'TBLA3'; ...
        'TBLA4'; ...
        'TBRC4'; ...
        'TBRC5'; ...
        'TBRB4'; ...
        'GC8'; ...
        'HL6'; ...
        'TLRA2'; ...
        'TBRA2'; ...
        'TBLB1'; ...
    };

     downsample_rate = 4;
     
     %records_to_process = 1:147;
     records_to_process = 1:9999999;
end  

%% Configuration for pat_FR_1125
if  strcmp(patient_id, 'pat_FR_1125')

    channels = { ...
        'HR12'; ...
        'HR10'; ...
        'HR11'; ...
        'HR13'; ...
        'HR14'; ...
        'HR9'; ...
        'TBA1'; ...
        'TBA2'; ...
        'TBA3'; ...
        'TBA4'; ...
        'TBB1'; ...
        'TBB2'; ...
        'TBB3'; ...
        'TBB4'; ...
        'TBC1'; ...
        'TBC2'; ...             
    };

     downsample_rate = 4;
     
     records_to_process = 81:199;
     %records_to_process = 123:130;
     %records_to_process = 6;
end  



%% Configuration for pat_FR_442
if  strcmp(patient_id, 'pat_FR_442')

    channels = { ...
        'HRA4'; ...
        'HRA5'; ...
        'TBA1'; ...
        'TBA2'; ...
        'TBA3'; ...
        'FRA1'; ...
        'FRA2'; ...
        'FRA3'; ...
        'FRA4'; ...
        'FRA5'; ...
        'FRA6'; ...
        'FRB1'; ...
        'FRB2'; ...
        'FRB3'; ...
        'FRB4'; ...
        'FRB5'; ...         
    };

     downsample_rate = 4;
     
     %records_to_process = 63:156;
     records_to_process = 69:91;
end  



end