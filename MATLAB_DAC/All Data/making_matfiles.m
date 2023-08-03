tic 

clc
clearvars

pi_co2 = 40;

price_data = readtable("WNSTRGN1_7_N001.csv");
emission_data = readtable("CAISO_full_filtered_emission.csv");
RTP = zeros(288,365);
for day = 1:365
    for row = 1:288
        %unadjusted
        RTP(row, day) = price_data{288*(day-1) + row, 4};
        
        %adjusted
        %RTP(row, day) = price_data{288*(day-1) + row, 4}+ pi_co2*emission_data{288*(day-1) + row +1, 11};
        
        toc
    end
end

%To save .mat file, right-click on the item "RTP" in the workspace