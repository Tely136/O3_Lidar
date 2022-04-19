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
n=floor(TimeLen/T);% new length of the time dimension
diff=mod(TimeLen,T);
DataTimeAve=zeros(lenHeight,n);
TimeIndArr=zeros(n,1);
for i=1:n% every 16 points average and store in a new matrix
    ind=(i-1)*T;
    DataTimeAve(:,i)=min(DataMat(:,ind+1:ind+T),[],2,'omitnan');
    TimeIndArr(i)=uint32(ind+floor((1+T)/2));
    if i==n
    DataTimeAve(:,i)=min(DataMat(:,ind+1:ind+T+diff),[],2,'omitnan');
    end
end