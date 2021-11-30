function [lpc_coeff, pitch_period, isvoiced, G,mn,delta,delay] = encoder_delay(arr, Fs,frame_length)

n= frame_length - 1;  

%Looping over frame to Calculate some important parameters
for frame_start=1 : frame_length : (length(arr) - frame_length)
    
    % Preemphasis to boost higher frequency of current frame
    y = filter([1 -.9378], 1, arr(frame_start:frame_start+n));
    
    % Magnitude Sum Calculation  (Formula used from theory)
    [first,second] = butter(9,.33,'low');  
    magnitude_sum(frame_start:(frame_start + n))=sum(abs(filter(first,second,y)));
    
    % Zero Crossing Calculation  (Formula used from theory)
    count=0;
    for n=1:length(y)-1
        count=count + (1./2) .* abs(sign(y(n+1))-sign(y(n)));
    end
    zero_crossing(frame_start:(frame_start + n)) =count;
    
    % Pitch Period Calculation   (Formula used from theory)
    femate_period = round (Fs .* 1e-3);
    male_period = round (Fs .* 25e-3);
    auto_correlation=xcorr(y);
    [a,bb]=max(auto_correlation); 
    pitch_per_range = auto_correlation ( bb + femate_period :bb + male_period );
    [a,bb]= max(pitch_per_range);
    pitch_period(frame_start:(frame_start + n)) = bb + femate_period;

end

%Making the boolean array on the basis of wheather the value is more than
%threshold value of the respective parameter
voiced_magnitude_sum_function =  check_threshhold(magnitude_sum,0.67,0);   
voiced_zero_crossing = check_threshhold(zero_crossing,1.5,1);
voiced_pitch_period =check_threshhold(pitch_period,0.5,0);

%Making the voiced/unvoiced array (1 when all 3 of the above are 1)
isvoiced=voiced_magnitude_sum_function.*voiced_pitch_period.* voiced_zero_crossing;

%Looping over to calculate gain and coefficients
for frame_start=1 : frame_length : (length(arr) - frame_length)
    
    % Preemphasis to boost higher frequency of current frame
    y = filter([1 -.9378], 1, arr(frame_start:frame_start+n));

    %Finding the Coefficient using DURBIN's recursion METHOD (10th order);
    buffer = dsp.AsyncBuffer(2*frame_length);
    write(buffer,y);
    hamming_window=hamming(2*frame_length);
    after_mult=hamming_window.*read(buffer,2*frame_length, frame_length);
    after_corr = xcorr(after_mult, 10, 'biased');
    after_corr = after_corr(11:end);
    [coefficient, ~, unused] = levinson(after_corr); 
    lpc_coeff(frame_start: (frame_start + length(coefficient) - 1))=coefficient; 
    
    %Error Calculation
    prediced_value= filter([0 -coefficient(2:end)],1,y);
    error_calculated= y -prediced_value;     

    %GAIN calculation   (Different formula for unvoiced and voiced frame from theory);
    if (~isvoiced(frame_start))
        G(frame_start)= sqrt( sum(error_calculated (1:length(error_calculated)) .^2) ./ length(error_calculated) );
    else
        G(frame_start)= sqrt( pitch_period(frame_start) .* sum( error_calculated (1:( floor( length(error_calculated)./pitch_period(frame_start) ) .* pitch_period(frame_start) )) .^2 ) ./ ( floor( length(error_calculated)./pitch_period(frame_start) ) .* pitch_period(frame_start) ) );
    end
end

% quantization in order to send the data in the form of bits
[lpc_coeff,mn(1),delta(1)]=quant(lpc_coeff,16);
[pitch_period,mn(2),delta(2)]=quant(pitch_period,16);
[isvoiced,mn(3),delta(3)]=quant(isvoiced,16);
[G,mn(4),delta(4)]=quant(G,16);

%Plotting of some important plots

