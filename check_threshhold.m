function [bool_output] = check_threshhold(arr,mult_factor,ismore)
 if(ismore==0) % Check if the sign should be more than or less than
    bool_output=arr > (( (sum(arr)./length(arr))-min(arr)) .*(mult_factor))+min(arr);  
 else
     bool_output=arr < (( (sum(arr)./length(arr))-min(arr)) .*(mult_factor))+min(arr);  
 end
