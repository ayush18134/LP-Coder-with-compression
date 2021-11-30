function output = decoder (mn,delta,lpc_coeff, pitch_period, isvoiced, G,frame_length);

%Converting the bitstream to data
lpc_coeff=bin2dec(lpc_coeff)';
pitch_period=bin2dec(pitch_period)';
isvoiced=bin2dec(isvoiced)';
G=bin2dec(G)';

% Dequantizing and getting the original values of coefficients, period,
% gain and isvoiced
lpc_coeff=mn(1)+lpc_coeff*delta(1);
pitch_period=floor(mn(2)+pitch_period*delta(2))+1;
isvoiced=mn(3)+isvoiced*delta(3);
G=mn(4)+G*delta(4);

n=length(lpc_coeff);

%Looping over the frame for decoding process
for frame_start=1:frame_length:n   
    %Checking if the frame is voiced or not 
    if (isvoiced(frame_start) == 1) % If voiced u will be an impulse train
        impulse_train=zeros(frame_length,1);
        for sample_num=1:frame_length
            if (rem(sample_num,pitch_period(frame_start))==0) %Making the array 1 at starting of each period
                impulse_train(sample_num) = 1;
            end
        end
        u=impulse_train;
    else
        awgn_noice= randn(1, frame_length); %If the frame is unvoiced u=AWGN noise
        u=awgn_noice;
    end
    % Calcultion of the reconstructed signal
    % Reconstructed signal= S(z)=H(z) *G *U(z)
    output(frame_start:frame_start+frame_length-1) = filter(1, [1 lpc_coeff((frame_start+1):(frame_start+1+9))], u).* G(frame_start);
end

%Plotting impulse train and AWGN noise
%figure;
%subplot(2,1,1), 
%plot(impulse_train ); 
%title("Impulse Train"); 
%subplot(2,1,2), 
%plot(awgn_noice); 
%title("Gaussian Noise");