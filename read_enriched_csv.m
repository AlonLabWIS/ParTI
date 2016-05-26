function [discreteTitles,discreteEnData] =  read_enriched_csv(fileName,delimiter)
  fid = fopen(fileName,'r');   
  discreteEn = cell(3000,1);     
                              
  line_ind = 1;               
  nextLine = fgetl(fid);       
  while ~isequal(nextLine,-1)        
    discreteEn{line_ind} = nextLine;  
    line_ind = line_ind+1;          
    nextLine = fgetl(fid);            
  end
  fclose(fid);                
  discreteEn = discreteEn(1:line_ind-1);  
  for i = 1:line_ind-1              
    data = textscan(discreteEn{i},'%s',...  
                        'Delimiter',delimiter);
    data = data{1};              
    if strcmp(discreteEn{i}(end),delimiter) 
      data{end+1} = '';                    
    end
    discreteEn(i,1:numel(data)) = data; 
  end
discreteTitles = discreteEn(1,:);
discreteEnData = discreteEn(2:end,:);
end
