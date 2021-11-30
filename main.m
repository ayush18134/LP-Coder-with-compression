%Ayush Bhardwaj
%2018134
%SAP Assignment 4
clear all;
clc;


%Reading the 8KHz input file
input = 'original_signal.wav';
[arr, Fs] =audioread(input); 

%Setting up the Frame size and length
frame_size =30e-3;    
frame_length = round(Fs .* frame_size);   

% Encoding the stream
[lpc_coeff, pitch_period, isvoiced, G ,mn,delta,processing_delay] = encoder(arr, Fs,frame_length);  

% Decoding the received signal
output = decoder(mn,delta,lpc_coeff, pitch_period, isvoiced, G,frame_length);


%RESULTS,
disp("Processing Delay for the encoder per frame");
disp(processing_delay);

%Save the file
audiowrite("reconstruced_signal.wav",output,Fs)
disp("The reconstructed signal audio file has been saved");

figure;
subplot(2,1,1), 
plot(arr); 
xlabel("Time in seconds");
ylabel("Amplitude");
title("Original Signal"); 
subplot(2,1,2), 
plot(output); 
xlabel("Time in seconds");
ylabel("Amplitude");
title("Reconstructed Signal");