function [scores, labels] = get_sepsis_score(input_file)
% read the data
load('method.mat');
Ta = ReadChallengeData(input_file);
load weights.mat
load hmmsimple.mat
D = Ta;
% discard the non significant features
if ~strcmp(method,'kmeans')
    D (:,indNonSig) = [];
    w(indNonSig) =[];
end

% processing Nans of the Data
if strcmp(method,'interp')
    D = DataInterpolationfcn(D, 1:size(D,2), 1);
else
    D(isnan(D))=0;
end
if ~strcmp(method,'kmeans')
    D = D./repmat(datanorm,[size(D,1),1]);
    D2 = D*w';
%     D2 = round(D2*N)./N;
%     seq=D2*N;
%     if min(seq)<=0
%         seq     = seq - min(seq)+1;
%     end
    D2 = (D2-min(D2))./( max(D2)-min(D2));
    seq =round(D2*N+1);
    
else
    for i=1:size(D,1)
        seq(i)=findCluster(centroids,D);
    end
    seq=seq';
end

PSTATES = hmmdecode(seq',trans_est,emis_est);
scores=PSTATES(2,:)';
% try
    labels = hmmviterbi(seq',trans_est,emis_est);
    
% catch e
% %     disp(e.message);
% end 
labels=labels-1;
labels = labels';
end

function [values, column_names] = ReadChallengeData(filename)
  f = fopen(filename, 'rt');
  try
    l = fgetl(f);
    column_names = strsplit(l, '|');
    values = dlmread(filename, '|', 1, 0);
  catch ex
    fclose(f);
    rethrow(ex);
  end
  fclose(f);

  %% ignore SepsisLabel column if present
  if strcmp(column_names(end), 'SepsisLabel')
    column_names = column_names(1:end-1);
    values = values(:,1:end-1);
  end
end
