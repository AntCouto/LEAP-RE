
clc; clear all;  tic;

%% Read RDA 461.0 database

DataFolder = './Algeria_2019';
% DataFolder = './Rwanda_Quenia_2020';
Outfolder = './FilesOutput';
EstFile = 'Algeria';
%% Processing
home = pwd;
addpath(home);
j=1;

Folderlist=dir('GDAS.ADPSFC.20*.txt'); % list all folders

cd(DataFolder)
[FolderLista,~]=size(Folderlist);
f = 1
% for f = 1: FolderLista
%     cd(Folderlist(f).name)
    
%     display(['Reading folder ', num2str(f), ' / of ', num2str(FolderLista)])
    Filesist=dir('gdas.ADPSFC.20*.txt'); % list all folders
    [lines,~]=size(Filesist);
     
    
    for k=1:lines % Cycle for read each file
     display(['Reading file ', num2str(k), ' / of ', num2str(lines)])
     
        openfile=fopen(Filesist(k).name); % open file
        line1=fgetl(openfile); % read 1st line of the file
        line2=fgetl(openfile); % read 2nd line of the file
        line3=fgetl(openfile); % read 3th line of the file
        
        while ~feof(openfile) % read each line until reach eof
            
            nextline=fgetl(openfile); % next line
            splits=strsplit(nextline,' '); % divide a linha consoante os espaços
            split=cell(splits);
            
            if str2double(split(8)) >0 & str2double(split(8))<360 & str2double(split(9)) >0 & str2double(split(9))<50
                
                
                rectype_IN(j)=split(1); % Id. of record type
                obstype_IN(j)=split(2); % Id. of observation type (e.g. SYNOP or METAR)
                time_IN(j)=datenum(split(3),'yyyymmddHHMM'); % convert to MATLAB date
                
                % Save date in a specific format:
                data_gravar1_IN(j) = str2num(datestr(time_IN(j),'yyyymmdd'));
                time_gravar1_IN(j) = str2num( datestr(time_IN(j),'HHMM')); clear time_gravar
                
                stationID_IN(j)=split(4);              % ID of the station
                latitude_IN(j)=str2double(split(5));  % Latitude
                longitude_IN(j)=str2double(split(6)); % Longitude
                elevation_IN(j)=str2double(split(7)); % Elevation
                
                winddirection_IN(j)=str2double(split(8)); % Wind direction
                windspeed_IN(j)=str2double(split(9)); %Wind speed
                
                j=j+1;
            end
        end
        closefile=fclose(openfile); % fecha o ficheiro
        
    end
    cd ..
% end
cd(Outfolder)

%% Ascendly sort the data :
[Time IdxOrd ] = sort(time_IN,'ascend');
 
whois = whos('*_IN');
for i = 1:length(whois);
    AA = eval(whois(i).name);
    BB = AA(IdxOrd);
     
    varname = strcat(whois(i).name);
    varnameF = varname(1:end-3);
    eval([varnameF '= BB ;']);
    clear(whois(i).name)
    clear AA BB varname
    
end

%%
clear f Filesist Folder* lines j whois
clear list  openfile line1 line2 line3 closefile k splits split nextline

%% Conversion of Estation ID to number
% More fast to work with number and some stations show letters

stationID_num = zeros(length(stationID),1);
for i = 1 : length(stationID)
    aux1 = str2num(cell2mat(stationID(i)));
    
    obstype_num (i,1) = sum(uint16(cell2mat(obstype(i))));
    if length(aux1) >= 1
        stationID_num (i,1) = aux1;
    else
        
        aux2 =  uint16(cell2mat(stationID(i)));
        if length(aux2) == 4
            stationID_num (i,1)  = 10*double(aux2(1))+1.3130*double(aux2(2))+0.71*double(aux2(3))+0.2333*double(aux2(4));
        else
            stationID_num (i,1)  = 10*double(aux2(1))+1.3130*double(aux2(2))+0.71*double(aux2(3))+0.2333*double(aux2(4))+ 1.0441*double(aux2(5));
        end
    end
    
end

% Contabilizar o número de registos por estação antes do controlo de
% qualidade
C=unique(stationID_num); % o vector C guarda o ID das estações existentes
H{1,1}='Station ID'; H{1,2}='Nr.observations';
H{1,3}='Latitude';  H{1,4}='Longitude';

for u=1:length(C)    % sem repetições
    
    pos = find ( C(u) == stationID_num);
    
    
    H(u+1,1)=stationID(pos(1)); % Na primeira coluna guarda o ID de cada estação
    H(u+1,2)=num2cell(length(pos)); % Na segunda coluna guarda o nº de ocorrências
    H(u+1,3) = num2cell(latitude(pos(1)));
    H(u+1,4) = num2cell(longitude(pos(1)));
    
    clear pos
    
end

clear count u C TF aux1 aux2

%% Criação dos ficheiros filtrados para a base de dados final
[CC , ia, ic] =unique(stationID_num); 
[CC1]= unique(obstype_num);

% Nome_Extenso=
ref_1 = sum(uint16(['SYNOP';'METAR';])');
ref_2 =  sum(uint16(['SYNOPM'])');
ref_3 = sum(uint16(['EFGHIJKL'])');
ref = [ref_1  ref_2 ref_3]; clear ref_*
Ref_Year = mode(year(Time));
Aux_Year = datenum(Ref_Year,1,1,0,0,0):datenum(0,0,0,1,0,0):datenum(Ref_Year+1,1,1,0,0,0);
%%
% cd (Outfolder)
for uu=1:length(CC)    
    clear media_ws30 media_ws1 media_wd30 media_wd1 aux30 aux1
    
    K = [];
    display(['Writting file ', num2str(uu) ,' of ', num2str(length(CC))])
    if (obstype_num( ia(uu)) == ref(1) | obstype_num( ia(uu)) == ref(2) | obstype_num( ia(uu)) == ref(1))
        
        aux1 = find ( CC(uu) == stationID_num);
        
        if  (obstype_num(aux1(1)) == ref(1) | obstype_num(aux1(1)) == ref(end))% Se for SYNOP retorna 1
            ID=stationID{aux1(1)};
            ID2=stationID{aux1(1)};
            % Header of the file
            string='_SYNOP_WIND';
            filename=strcat(ID2,string); % Nome do ficheiro
            K{1,1}='RecType'; K{1,2}='ObsType'; K{1,3}='Matlab Time';
            K{1,4}='Latitude(º)'; K{1,5}='Longitude(º)';
            K{1,6}='Elevation(m)'; K{1,7}='WindDirection(º)';
            K{1,8}='WindSpeed(m/s)';

            K(2:length(aux1)+1,1)=rectype(aux1); % Preenchimento do ficheiro
            K(2:length(aux1)+1,2)=obstype(aux1);
            K(2:length(aux1)+1,3)=num2cell(time(aux1));
            K(2:length(aux1)+1,4)=num2cell(latitude(aux1));
            K(2:length(aux1)+1,5)=num2cell(longitude(aux1));
            K(2:length(aux1)+1,6)=num2cell(elevation(aux1));
            K(2:length(aux1)+1,7)=num2cell(winddirection(aux1));
            K(2:length(aux1)+1,8)=num2cell(windspeed(aux1));
            
            
            [Ti Ta Tb] = intersect(Aux_Year,time(aux1));
            EstFin{uu,1}=ID;
            EstFin{uu,2}=latitude(aux1(1));
            EstFin{uu,3}=longitude(aux1(1));
            EstFin{uu,4}=length(aux1);
            EstFin{uu,5}=100*(length(Ti)/length(Aux_Year));
            
            nomeID= (CC(uu));
            
            nomeID2=  str2num(cell2mat(stationID(ia(uu))));
            data_ini=datevec(time(1)+1); % soma-se 1 dia para garantir que se está no mês correto e evitar que possa estar no fim do mês anterior
            ano_i=num2str(data_ini(1));
            if data_ini(2)<10
                mes_i=num2str(data_ini(2),'%02.0f');
            else
                mes_i=num2str(data_ini(2));
            end
            
            if nomeID == nomeID2
                
                nomeID = num2str(nomeID);
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], K,'DataBase','A1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'ID_Station:'},'DataBase','k1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],  {nomeID},'DataBase','l1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Hour'},'DataBase','i1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Date'},'DataBase','j1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],[time_gravar1(aux1)' data_gravar1(aux1)' ],'DataBase','i2');
                
            else
                nomeID = cell2mat(stationID(ia(uu)));
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], K,'DataBase','A1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'ID_Station:'},'DataBase','k1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],  {nomeID},'DataBase','l1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Hour'},'DataBase','i1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Date'},'DataBase','j1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],[time_gravar1(aux1)' data_gravar1(aux1)' ],'DataBase','i2');
            end
            clear K
            clear data_ini ano_i mes_i
        end
    
        if  obstype_num(aux1(1)) == ref(3) % Se for SYNOP retorna 1
            
            ID=stationID{aux1(1)};
            ID2=stationID(ia(uu));
            % Header of the file
            string='_SYNOPM_WIND';
            filename=strcat(ID,string); % Nome do ficheiro
            K{1,1}='RecType'; K{1,2}='ObsType'; K{1,3}='Matlab Time';
            K{1,4}='Latitude(º)'; K{1,5}='Longitude(º)';
            K{1,6}='Elevation(m)'; K{1,7}='WindDirection(º)';
            K{1,8}='WindSpeed(m/s)';
            
            K(2:length(aux1)+1,1)=rectype(aux1); % Preenchimento do ficheiro
            K(2:length(aux1)+1,2)=obstype(aux1);
            K(2:length(aux1)+1,3)=num2cell(time(aux1));
            K(2:length(aux1)+1,4)=num2cell(latitude(aux1));
            K(2:length(aux1)+1,5)=num2cell(longitude(aux1));
            K(2:length(aux1)+1,6)=num2cell(elevation(aux1));
            K(2:length(aux1)+1,7)=num2cell(winddirection(aux1));
            K(2:length(aux1)+1,8)=num2cell(windspeed(aux1));
            
            [Ti, ~, ~] = intersect(Aux_Year,time(aux1));
            EstFin{uu,1}=ID;
            EstFin{uu,2}=latitude(aux1(1));
            EstFin{uu,3}=longitude(aux1(1));
            EstFin{uu,4}=length(aux1);
            EstFin{uu,5}=100*(length(Ti)/length(Aux_Year));
            
            nomeID= (CC(uu));
            
            nomeID2=  str2num(cell2mat(stationID(ia(uu))));
            data_ini=datevec(time(1)+1); % soma-se 1 dia para garantir que se está no mês correto e evitar que possa estar no fim do mês anterior
            ano_i=num2str(data_ini(1));
            if data_ini(2)<10
                mes_i=num2str(data_ini(2),'%02.0f');
            else
                mes_i=num2str(data_ini(2));
            end
            
            if nomeID == nomeID2
                nomeID = num2str(nomeID);
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], K,'DataBase','A1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'ID_Station:'},'DataBase','k1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],  {nomeID},'DataBase','l1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Hour'},'DataBase','i1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Date'},'DataBase','j1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],[time_gravar1(aux1)' data_gravar1(aux1)' ],'DataBase','i2');
                
                
            else
                nomeID = cell2mat(stationID(ia(uu)));
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], K,'DataBase','A1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'ID_Station:'},'DataBase','k1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],  {nomeID},'DataBase','l1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Hour'},'DataBase','i1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Date'},'DataBase','j1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],[time_gravar1(aux1)' data_gravar1(aux1)' ],'DataBase','i2');
            end
            clear K
            clear data_ini ano_i mes_i
        end
        
        
        % Para os ficheiros com observações METAR
        if obstype_num(aux1(1)) == ref(2)          % Se for METAR retorna 0
            ID=stationID{aux1(1)};
            ID2=stationID(ia(uu));
            string='_METAR_WIND';
            filename=strcat(ID,string);
            K{1,1}='RecType'; K{1,2}='ObsType'; K{1,3}='Matlab Time'; ...
                K{1,4}='Latitude(º)'; K{1,5}='Longitude(º)'; ...
                K{1,6}='Elevation(m)'; K{1,7}='WindDirection(º)'; ...
                K{1,8}='WindSpeed(m/s)';
            
            
            K(2:length(aux1)+1,1)=rectype(aux1); % Preenchimento do ficheiro
            K(2:length(aux1)+1,2)=obstype(aux1);
            K(2:length(aux1)+1,3)=num2cell(time(aux1));
            K(2:length(aux1)+1,4)=num2cell(latitude(aux1));
            K(2:length(aux1)+1,5)=num2cell(longitude(aux1));
            K(2:length(aux1)+1,6)=num2cell(elevation(aux1));
            K(2:length(aux1)+1,7)=num2cell(winddirection(aux1));
            K(2:length(aux1)+1,8)=num2cell(windspeed(aux1));
            
            
            [Ti, ~, ~] = intersect(Aux_Year,time(aux1));
            EstFin{uu,1}=ID;
            EstFin{uu,2}=latitude(aux1(1));
            EstFin{uu,3}=longitude(aux1(1));
            EstFin{uu,4}=length(aux1);
            EstFin{uu,5}=100*(length(Ti)/length(Aux_Year));
            
            nomeID= (CC(uu));
            
            nomeID2=  str2num(cell2mat(stationID(ia(uu))));
            data_ini=datevec(time(1)+1); % soma-se 1 dia para garantir que se está no mês correto e evitar que possa estar no fim do mês anterior
            ano_i=num2str(data_ini(1));
            if data_ini(2)<10
                mes_i=num2str(data_ini(2),'%02.0f');
            else
                mes_i=num2str(data_ini(2));
            end
            
            if nomeID == nomeID2
                
                nomeID = num2str(nomeID);
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], K,'DataBase','A1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'ID_Station:'},'DataBase','k1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],  {nomeID},'DataBase','l1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Hour'},'DataBase','i1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Date'},'DataBase','j1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],[time_gravar1(aux1)' data_gravar1(aux1)' ],'DataBase','i2');
                
            else
                
                nomeID = cell2mat(stationID(ia(uu)));
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], K,'DataBase','A1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'ID_Station'},'DataBase','k1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],  {nomeID},'DataBase','l1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Hour'},'DataBase','i1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'], {'Date'},'DataBase','j1');
                xlswrite(['DataBase_ID_',nomeID,'_',mes_i,'_',ano_i,'.xlsx'],[time_gravar1(aux1)' data_gravar1(aux1)' ],'DataBase','i2');
                
            end
            
            clear K
            clear data_ini ano_i mes_i
        end
        
    end
end

%% Save some statistics
xlswrite([EstFile,'_',num2str(Ref_Year),'.xlsx'], [{'Station'},{'Lat.'},{'Lon.'},{'TotalRec.'},{'1hrAveAva..'}],'Stat.','A1');
xlswrite([EstFile,'_',num2str(Ref_Year),'.xlsx'], EstFin,'Stat.','A2');


%%
% Representa as observações na área escolhida incluindo a linha de costa
fig2=figure('name','Disponibilidade Espacial de Estações');
load coast
plot(long,lat,'k')

hold on
scatter(cell2mat(H(2:end,4)),cell2mat(H(2:end,3)),30,cell2mat(H(2:end,2)),'Filled');
xlim([min(cell2mat(H(2:end,4)))-1 max(cell2mat(H(2:end,4)))+1])
ylim([min(cell2mat(H(2:end,3)))-1 max(cell2mat(H(2:end,3)))+1])
title('\bfNumber of records')
xlabel('Longitude (º)'); ylabel('Latitude (º)')
colorbar


toc