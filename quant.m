function [out,mn,delta] = quant(arr,bits)
mx=max(arr);
mn=min(arr);
% Definig width of each level
delta=(mx-mn)/(2^bits);
out=zeros(size(arr));
for i=1:length(out)
    %Using standard Formula for the process of quantization
    out(i)=floor((arr(i)-mn)/delta);
end
end
