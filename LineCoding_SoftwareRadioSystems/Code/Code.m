clc
clear
%Amplitude the encoded bits
A=4;
%Number of Realizations
numberofRealizations=500;
%Total Number of Samples
samplesPerRealization=700;
%Number of Generated Bits
numberofBits=100;
%Number of Samples per Bit
samplesPerBit=7;
%Bit Duration
Tb=0.07;
%Sample Duration
Ts=0.01;
%Sampling Frequency
fs =1/Ts;

%Generating random bitstream of 101 bits for 500 Realizations
Data=randi([0 1],numberofRealizations,numberofBits+1);
%Generating different initial time delays for each Realization
T_delay=randi([0 6],numberofRealizations,1);
for i=1:3
  if i==1
  %For UniPolar 
  Tx_UniPolar=(Data)*A;
  %Tx(500x101) Matrix where the each element in Tx is repeated 7 times
  Tx2_UniPolar=repelem(Tx_UniPolar,samplesPerBit,1);
  %Reshaping the Tx Matrix into a single Column
  %each 707 Row of Tx_out represent a Realization
  Tx_out_UniPolar=reshape(Tx2_UniPolar,size(Tx2_UniPolar,1)*size(Tx2_UniPolar,2),1);
  %each Column will represent a Realization 
  Tx_out_UniPolar=reshape(Tx_out_UniPolar,(numberofBits+1)*samplesPerBit,numberofRealizations);
  %Matrix Transpose each column into a row , each row into a column
  %Now each Row Represent a Realization     
  Tx_out_UniPolar=Tx_out_UniPolar';
  %Adding Random Initial Time Delay
  Tx_out_new_UniPolar=time_delay(Tx_out_UniPolar,T_delay);
  %Computing Statistical Mean
  statistical_mean_UniPolar=stat_mean(Tx_out_new_UniPolar);
  %Computing Statistical AutoCorrelation          
  statistical_autocorrelation_UniPolar=statistical_autocorr(Tx_out_new_UniPolar);
  %Computing Time Mean
  time_mean_UniPolar=time_mean_array(Tx_out_new_UniPolar);
  %Computing Time Auto Correlation
  time_autocorrelation_UniPolar=time_autocorrelation_array(Tx_out_new_UniPolar);
    elseif i==2
  %For Polar NRZ
  Tx_Polar_NRZ=((2*Data)-1)*A;
  %Tx(500x101) Matrix where the each element in Tx is repeated 7 times
  Tx2_Polar_NRZ=repelem(Tx_Polar_NRZ,samplesPerBit,1);
  %Reshaping the Tx Matrix into a single Column
  %each 707 Row of Tx_out represent a Realization
  Tx_out_Polar_NRZ=reshape(Tx2_Polar_NRZ,size(Tx2_Polar_NRZ,1)*size(Tx2_Polar_NRZ,2),1);
  %each Column will represent a Realization 
  Tx_out_Polar_NRZ=reshape(Tx_out_Polar_NRZ,(numberofBits+1)*samplesPerBit,numberofRealizations);
  %Matrix Transpose each column into a row , each row into a column
      %Now each Row Represent a Realization 
      Tx_out_Polar_NRZ=Tx_out_Polar_NRZ';
      %Adding Random Initial Time Delay
      Tx_out_new_Polar_NRZ=time_delay(Tx_out_Polar_NRZ,T_delay);
      %Computing Statistical Mean
      statistical_mean_Polar_NRZ=stat_mean(Tx_out_new_Polar_NRZ);
      %Computing Statistical AutoCorrelation
      statistical_autocorrelation_Polar_NRZ=statistical_autocorr(Tx_out_new_Polar_NRZ);
       %Computing Time Mean
        time_mean_Polar_NRZ=time_mean_array(Tx_out_new_Polar_NRZ);
        %Computing Time Auto Correlation
        time_autocorrelation_Polar_NRZ=time_autocorrelation_array(Tx_out_new_Polar_NRZ);
   else
      %For Polar RZ
      Tx_Polar_RZ=((2*Data)-1)*A;
      %Tx(500x101) Matrix where the each element in Tx is repeated 7 times
      Tx2_Polar_RZ=repelem(Tx_Polar_RZ,samplesPerBit,1);
      %Reshaping the Tx Matrix into a single Column
      %each 707 Row of Tx_out represent a Realization
      Tx_out_Polar_RZ=reshape(Tx2_Polar_RZ,size(Tx2_Polar_RZ,1)*size(Tx2_Polar_RZ,2),1);
      %For Polar RZ first 4 samples of each bit have value A or -A else 0
      for l=5:samplesPerBit:size(Tx_out_Polar_RZ,1)
          for j=l:1:l+2
             Tx_out_Polar_RZ(j)=0;
          end
      end     
        %each Column will represent a Realization 
        Tx_out_Polar_RZ=reshape(Tx_out_Polar_RZ,(numberofBits+1)*samplesPerBit,numberofRealizations);
        %Matrix Transpose each column into a row , each row into a column
        %Now each Row Represent a Realization 
        Tx_out_Polar_RZ=Tx_out_Polar_RZ';
        %Adding Random Initial Time Delay
        Tx_out_new_Polar_RZ=time_delay(Tx_out_Polar_RZ,T_delay);
        %Computing Statistical Mean
        statistical_mean_Polar_RZ=stat_mean(Tx_out_new_Polar_RZ);
        %Computing Statistical AutoCorrelation
        statistical_autocorrelation_Polar_RZ=statistical_autocorr(Tx_out_new_Polar_RZ);
        %Computing Time Mean
        time_mean_Polar_RZ=time_mean_array(Tx_out_new_Polar_RZ);
         %Computing Time Auto Correlation
        time_autocorrelation_Polar_RZ=time_autocorrelation_array(Tx_out_new_Polar_RZ);
    end
end    
  %Plotting Line Codes
    figure();
    subplot(3,1,1);
    plot_line_code(Tx_out_new_UniPolar)
    title("UniPolar Line Code");
    subplot(3,1,2);
    plot_line_code(Tx_out_new_Polar_NRZ)
    title("Polar NRZ Line Code");
    subplot(3,1,3);
    plot_line_code(Tx_out_new_Polar_RZ)
    title("Polar RZ Line Code");
   %Plotting Statistical mean
    figure();
    subplot(3,1,1);
    plot_statistical_mean(statistical_mean_UniPolar,A);
    title("UniPolar Statistical Mean");
    subplot(3,1,2);
    plot_statistical_mean(statistical_mean_Polar_NRZ,A);
    title("Polar NRZ Statistical Mean");
    subplot(3,1,3);
    plot_statistical_mean(statistical_mean_Polar_RZ,A);
    title("Polar RZ Statistical Mean");
    %Plotting Statistical Autocorrelation
    figure();
    subplot(3,1,1);
    plot_statistical_autocorr(statistical_autocorrelation_UniPolar);
    title("Statistical AutoCorrelation Function of UniPolar");
    ylim([0,10]);
    subplot(3,1,2);
    plot_statistical_autocorr(statistical_autocorrelation_Polar_NRZ);
    title("Statistical AutoCorrelation Function of Polar NRZ");
    ylim([-2,20]);
    subplot(3,1,3);
    plot_statistical_autocorr(statistical_autocorrelation_Polar_RZ);
    title("Statistical AutoCorrelation Function of Polar RZ");
    ylim([0,10]);
    %Plotting Time Mean
    figure();
    subplot(3,1,1);
    plot_time_mean(time_mean_UniPolar,A)
    title("Time Mean of UniPolar");
    subplot(3,1,2);
    plot_time_mean(time_mean_Polar_NRZ,A)
    title("Time Mean of Polar NRZ");
    subplot(3,1,3);
    plot_time_mean(time_mean_Polar_RZ,A)
    title("Time Mean of Polar RZ");
    %Plotting Time Autocorrelation
    figure();
    subplot(3,1,1);
    plot_time_autocorr(time_autocorrelation_UniPolar);
    title("Time Autocorrelation of UniPolar");
    ylim([0,10]);
    subplot(3,1,2);
    plot_time_autocorr(time_autocorrelation_Polar_NRZ);
    title("Time Autocorrelation of Polar NRZ");
    ylim([-2,20]);
    subplot(3,1,3);
    plot_time_autocorr(time_autocorrelation_Polar_RZ);
    title("Time Autocorrelation of Polar RZ");
    ylim([0,10]);
    %Plotting Power Spectral Density
    figure();
    subplot(3,1,1);
    plot_psd(statistical_autocorrelation_UniPolar);
    title("PSD of UniPolar");
    hold on
    plot_PSD_theoretical_UniPolar(A,Tb)
    ylim([0,5]);
    legend('Actual','Theoretical');
    subplot(3,1,2);
    plot_psd(statistical_autocorrelation_Polar_NRZ);
    title("PSD of Polar NRZ");
        hold on
plot_PSD_theoretical_polar_NRZ(A,Tb)
    ylim([0,0.3]);
        legend('Actual','Theoretical');
    subplot(3,1,3);
    plot_psd(statistical_autocorrelation_Polar_RZ);
    title("PSD of Polar RZ");
    hold on
    plot_PSD_theoretical_polar_RZ(A,Tb)
    ylim([0,0.3]);
    legend('Actual','Theoretical');
   
 function output_array=statistical_autocorr(Tx_out_new)
  [rows,cols]=size(Tx_out_new);
   output_array = zeros(1,cols);
   for col=1:cols
    sum=0;
        for row=1:rows
        sum=sum+Tx_out_new(row,1)*Tx_out_new(row,col);
        end
    mean=sum/rows;
    output_array(1,col)=mean;
   end
 end
 
 function output_array=stat_mean(Tx_out_new)
    [rows,cols]=size(Tx_out_new);
    output_array=zeros(1,cols);
  for col=1:cols
    sum=0;
        for row=1:rows
        sum=sum+Tx_out_new(row,col);
        end
  mean=sum/rows;
  output_array(1,col)=mean;
  end
 end
 function output_array=time_mean_array(Tx_out_new)
    [rows,cols]=size(Tx_out_new);
    output_array=zeros(1,rows);
    for row=1:rows
    sum=0;
        for col=1:cols
        sum=sum+Tx_out_new(row,col);
        end
    output_array(1,row)=sum/cols;
    end
 end
 function output_array=time_autocorrelation_array(Tx_out_new)
   [rows,cols]=size(Tx_out_new);
   output_array=zeros(1,rows);
   shift=1;
    for row=1:rows
    sum=0;
    shifted_wave=circshift(Tx_out_new(row,:),[0,shift-1]);
        for col=1:cols
            sum=sum+Tx_out_new(row,col)*shifted_wave(1,col);
        end
    output_array(1,row)=sum/cols;
    shift=shift+1;
    end
 end
 function output_array=time_delay(Tx_out,T_delay)
 samplesPerRealization=700;
 numberofRealizations=500;
 output_array=zeros(1,samplesPerRealization);
  %Adding Random Initial Time Delay for each Realization
    for count=1:numberofRealizations
        start_index=T_delay(count,1);
        end_index=(samplesPerRealization-1)+start_index;
  %Taking 700 samples for each Realization starting from the initial Delay
        sliced_row=Tx_out(count,start_index+1:end_index+1);
        output_array(count,:)=sliced_row;
    end
 end
 function plot_line_code(Tx_out_new)
 samplesPerBit=7;
 plot(Tx_out_new(1,1:samplesPerBit*7),'LineWidth',2)
 xlabel("time")
 ylabel("Magnitude")
 end
 function plot_statistical_mean(stat_mean,A)
 %Plot Statistical Mean
    samplesPerRealization=700;
    x = 1:1:samplesPerRealization;
    y = stat_mean(1,x);
    plot(x, y,'LineWidth', 2)
    ylim([-A, A])
    hold on
    xlabel("Samples")
    ylabel("Statistical Mean")
 end
 function plot_statistical_autocorr(stat_autocorr)
    %Plot Statistical AutoCorrelation
    samplesPerRealization=700;
    Ts=0.01;
    x=(-samplesPerRealization+1:1:samplesPerRealization)*Ts*1000;
    left_autocorrelation_mean=fliplr(stat_autocorr);
    full_autocorrelation_mean=horzcat(left_autocorrelation_mean,stat_autocorr);
    plot(x,full_autocorrelation_mean);
    hold on
    xlabel("tau[ms]")
    ylabel("Statistical AutoCorrelation")
 end
 function plot_time_mean(time_mean,A)
 %Plot Time Mean
   numberofRealizations=500;
   x=1:1:numberofRealizations;
   plot(x,time_mean,'LineWidth',2);
   ylim([-A,A]);
   hold on
   xlabel("Samples")
   ylabel("Time Mean")
 end
 function plot_time_autocorr(time_autocorrelation)
  %Plot Time AutoCorrelation 
    numberofRealizations=500;
    x=-numberofRealizations+1:1:numberofRealizations;
    left_time_autocorrelation=fliplr(time_autocorrelation);
    full_time_autocorrelation=horzcat(left_time_autocorrelation,time_autocorrelation);
    plot(x,full_time_autocorrelation);
    hold on
    xlabel("Samples")
    ylabel("Time Autocorrelation")
 end
 function plot_psd(stat_autocorr)
    Ts=0.01;
    fs = 1/Ts;
    samplesPerRealization=700;
    left_autocorrelation_mean=fliplr(stat_autocorr);
    full_autocorrelation_mean=horzcat(left_autocorrelation_mean,stat_autocorr);
   %Plot Power Spectral Density
   %Fourier Transform to obtain PSD
   normalized_psd=fft(full_autocorrelation_mean)/fs;
   conj_normalized_psd=conj(normalized_psd);
   squared_psd=conj_normalized_psd.*normalized_psd;
   square_root_psd=sqrt(squared_psd);
   psdshift=fftshift(square_root_psd);
   f = (-samplesPerRealization+1:1:samplesPerRealization)*fs/length(normalized_psd);
   plot(f,abs(psdshift),'LineWidth',2);
    hold on
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
 end
 function plot_PSD_theoretical_UniPolar(A,T)
 Ts=0.01;
 fs=1/Ts;
 f =(-699:1:700)*fs/700;
   % Calculate PSD components
    DC_component = (A^2 / 4) * (f == 0); % DC component
    sinc_component = (A^2 / (4 * T)) * (T * sinc(pi * f * T)).^2; % Sinc component

    % Total PSD
    S = DC_component + sinc_component;
    % Plot PSD
    plot(f, S,'LineWidth',2);
    
 end
 function plot_PSD_theoretical_polar_NRZ(A,T)
    Ts=0.01;
    fs=1/Ts;
    f =(-699:1:700)*fs/700;
    S=((A*A))*T*sinc(pi*f*T).^2;
    plot(f,S,'LineWidth',2);
 end
 function plot_PSD_theoretical_polar_RZ(A,T)
   Ts=0.01;
   fs =1/Ts;
    f =(-699:1:700)*fs/700;
    S=(((A*A))/4)*T*sinc(pi*f*T/2).^2;
    plot(f,S,'LineWidth',2);
 end