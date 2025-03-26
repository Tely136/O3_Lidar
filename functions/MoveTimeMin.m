%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling using linear interpolation 
% Input: 
%        DataMat double  2D matrix first dimension is height and second dimension is time
%        lenHeight double  Length of the height dimension 
%        TimeLen   double  Length of the time dimension
%        T         single  num of data need to be averaged
%
% Output: DataTimeAve double  Matrix after time averag and downsampling
%         TimeIndArr  uint32  Time index Array of the new matrix    
function [DataTimeAve,TimeIndArr]=MoveTimeMin(DataMat,lenHeight,TimeLen,T)
n=round(TimeLen/T);% new length of the time dimension
diff=n*T-TimeLen; 
DataTimeAve=zeros(lenHeight,n);
TimeIndArr=zeros(n,1);
if diff >0 % for example, if TimeLen =75, and T =10; then n round up to 8: diff= 80-75 = 5 >0
    for i=1:n% every n points average and store in a new matrix
        ind=(i-1)*T;
        if i<n % the first 7 groups are the average of 10 numbers 
            DataTimeAve(:,i)=min(DataMat(:,ind+1:ind+T),[],2,'omitnan');
            TimeIndArr(i)=uint32(ind+floor((1+T)/2));
        else % the last 8th is the average of the last 5 
            DataTimeAve(:,i)=min(DataMat(:,ind+1:ind+diff),[],2,'omitnan');
            TimeIndArr(i)=uint32(ind+floor((1+diff)/2));
        end
    end
else % otherwise if TimeLen = 72, n round down to 7, diff = 70-72 <=0
    for i=1:n% every n points average and store in a new matrix
        ind=(i-1)*T;
        DataTimeAve(:,i)=min(DataMat(:,ind+1:ind+T),[],2,'omitnan');
        TimeIndArr(i)=uint32(ind+floor((1+T)/2));
    end
end