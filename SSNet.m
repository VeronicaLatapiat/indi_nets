%**************Part 1:Input the information of samples and network information****
%**************sample information**************
%Example:TCGA-example cancer data (BRCA cancer datasets)
%expression_normal_fileName = 'Example_n.txt';

function ssn = SSNet(expression_normal_fileName)
expression_normal_fileName = '/shared/D7/vlatapiat/Data_simulada/Caso_A/TO-METHODS_DS-to-PCA_CasoA.txt';
%expression_normal_fileName = 'simulated_expression_data.txt';

%we output the sample-specific driver profiles by using different control
%methods
%   Input:
%         expression_fileName including expression_tumor_fileName and expression_normal_fileName)   
%         index:denotes we use which network construction method

%  Output:
%         The sample-specific driver profiles of PNC;
%         The column is the samples and the rows is the genes. The value ?1? denoted that the gene is driver genes; 
%************************part1:LOAD sample data and network data************************
%********************obtain the paired expression data******************

[normal,~,name_normal]=importdata(expression_normal_fileName);
gene_list=normal.textdata(2:end,1);Sample_name_normal=normal.textdata(1,2:end);normal_data=normal.data;
ref_data=normal_data;

 %*************** SSN *********************** 
% devuelve el número de filas y columnas cuando A es una matriz.
[n1,n2] = size(ref_data);
%disp(n1)
%disp(n2)

%disp(gene_list)
fileID = fopen('Gene_list_SSN.txt','w');
fprintf(fileID, '%s\n', gene_list{:}); %cell array of strings
fclose(fileID);


  for i=1:n2
   
    disp(i)
    sample_red=ref_data(:,i);
    [R0,P]=SSN(sample_red,ref_data);
    %   Output:
    %   adjacency_matrix:the network structure
    
    %disp(R0)
    %disp(P)
    
    [row,col] = size(R0);
      
      
        for j=1:col
            for k=1:row
                
                if abs(R0(j,k))>= 0.001
                %if abs(R0(j,k))>= 0.001
                %if abs(R0(j,k))>= 0.01 %tal punto de la matriz dado por j(col) y k(row)
                    R0(j,k)=0;% si es mayor o = a 0.01, se coloca un 0 en la red.
                    %disp(R0(j,k))
                elseif j == k
                    R0(j,k)=0;
                else % < 0.001
                    R0(j,k)=1;
                    %disp(R0(j,k))
                end
            end
            
        end
        
        formatSpec = "%s_%d";
        A3 = Sample_name_normal{i};
        A1 = i;
        str = sprintf(formatSpec,A3,A1);
        writematrix(R0,str,'Delimiter','tab')         
 	%dlmwrite(str, R0,'delimiter', '\t')	
  end
 
