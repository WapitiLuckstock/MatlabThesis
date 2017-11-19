function Dijlenmolens_GUI
%% UI Generation
f = figure('Visible','off','units','normalized','outerposition',[0,0,1,1],'Numbertitle','off');
hbutton1 = uicontrol('Style','pushbutton',...
    'String','Stooklijn','Position',[1000,600,180,25],...
    'Callback',@hbutton1_Callback);
hbutton2 = uicontrol('Style','pushbutton',...
    'String','Warmtevraag','Position',[1000,590,180,25],...
    'Callback',@hbutton2_Callback);
hbutton3 = uicontrol('Style','pushbutton',...
    'String','Elektriciteitsvraag','Position',[1000,580,180,25],...
    'Callback',@hbutton3_Callback);
hbutton4 = uicontrol('Style','pushbutton',...
    'String','Waterrad Model','Position',[1000,570,180,25],...
    'Callback',@hbutton4_Callback);
hbutton5 = uicontrol('Style','pushbutton',...
    'String','Inkomsten','Position',[1000,560,180,25],...
    'Callback',@hbutton5_Callback);
hbutton6 = uicontrol('Style','pushbutton',...
    'String','Variabele Kosten','Position',[1000,550,180,25],...
    'Callback',@hbutton6_Callback);
hbutton7 = uicontrol('Style','pushbutton',...
    'String','Piekduur','Position',[1000,540,180,25],...
    'Callback',@hbutton7_Callback);
hbutton8 = uicontrol('Style','pushbutton',...
    'String','Buffervat','Position',[1000,510,180,25],...
    'Callback',@hbutton8_Callback);
ha = axes('Units','pixels','Position',[50,100,900,800]);
f.Units = 'normalized';
ha.Units = 'normalized';
align([hbutton1,hbutton2,hbutton3,hbutton4,hbutton5,hbutton6,hbutton7,hbutton8],'Center','Fixed',7);
%% Data Generation
% excelfiles inladen
verbruik=xlsread('VerbruiksprofielDH2003','sheet1','B2:AW366');
Temp=xlsread('Temperaturen2003','blad1','E2:AZ366');
Radiatie=xlsread('2016Instralingsflux.xls','blad1','B2:AW366');
Warmtevraag = xlsread('36AppWarmtevraag');

%gegevens en aannames
Verwarmingskost=20572.99; % doorgegeven via mail
Aantalapp=36; % niet zeker
Prijsgas=0.03280; % exclusief btw https://eneco.be/nl/elektriciteit-en-gas/prijs-gas
Verliezenoudnet=0.4; % aanname
Correctievoordrukkemomentenaanvoer=5; % +5-5 op aanvoer
Minbuiten=0;Maxbuiten=16;Minaanvoer=55;Maxaanvoer=70; % buitentemp -> aanvoertemp
warmtepompinout=18;
waterradgemoutput=20; %kwh
waterradplusmin=0.40;
exportgas=0.035; % aanname
prijswaterrad=5.500*waterradgemoutput; % aanname
totalekost=200000; %aanname
Areazonnecollector=10; % variabele aanname
c_p = 4.18; %kJ/kgK warmtecoefficient van water
fillrate = 0.5; %Liter/s
V_tank = 2000; %Tank Volume in Liter
t_filled = V_tank/fillrate; %tijd tot de tank volledig gevuld is
t=0:1:t_filled;
water_volume_in_tank = fillrate*t; %ogenblikkelijke waterlevel
T_storage = 273+70; %Afgiftewarmte
T_discharged = 273+50; %Terugkeerwarmte
Warmtepompvermogen = 18*ones(365,24); %Productieprofiel van warmtepomp, deze gaat nog nauwkeuriger worden met waterradproductieprofiel
cteCOP = 3.6; %berekende sCOP
W_water = water_volume_in_tank*c_p*(T_storage-T_discharged)*(1/3600); %opgeslagen warmte in kWh
W_tank = max(W_water); %maximaal opgslaagbare warmte in de tank

%Tijdvector
firstDay = datetime(2016,1,1);
lastDay = datetime(2016,12,31);
Datevector = firstDay:caldays(1):lastDay;
time=ones(365,24);
for days = 1:365
    time(days,:) = 0:1:23;
end
time = hours(time);

% stooklijn
Aanvoertemp=[45 55 70 75];
COP=[5.2 4.06 2.85 2.32];
b=polyfit(Aanvoertemp,COP,1); %y=ax+b -> b(1) = a, b(2) = b
X=1:80;

% verband invoertemp dijle tov uitvoertemp met de COP ( gegevens van
% warmtepomp)
Aanvoertempdijle=[13.2 12.6 12 11.5 10.3 9.8 9.2 8.6];
Leavingtempdijle=7.5;
Deltatempdijle=Aanvoertempdijle-ones(size(Aanvoertempdijle))*Leavingtempdijle;


%berekening van aangenomen verlies
[m,n]=size(verbruik);
Aantalvakjes=m*n;
Voorbeeldverbruik=sum(sum(verbruik));
Vraagperapp=(Verwarmingskost/Aantalapp)/Prijsgas;
Correctietoevoegingwarmte=Verliezenoudnet*Voorbeeldverbruik/Aantalvakjes;
Pereenheidcorrectie=Vraagperapp/(Voorbeeldverbruik+Correctietoevoegingwarmte*Aantalvakjes);
% omzetting van uurlijk totaalverbruik naar werkelijk uurlijk verbruik zonder verlies
Echtverbruik=verbruik*Pereenheidcorrectie;

%correctiefactor op warmtevraag berekening
Gemiddeldewaardeperhalfuur=mean(Echtverbruik);
Maxgemtwaardeperhalfuur=max(Gemiddeldewaardeperhalfuur);
Mingemtwaardeperhalfuur=min(Gemiddeldewaardeperhalfuur);
[m2,n2]=size(Gemiddeldewaardeperhalfuur);
for i = 1:n2
    Gemiddeldewaardeperhalfuur(1,i)=(Gemiddeldewaardeperhalfuur(1,i)-Mingemtwaardeperhalfuur)/(Maxgemtwaardeperhalfuur-Mingemtwaardeperhalfuur)*10-5;
end

% buitentemp naar aanvoertemp
Aanvoertempvbuiten=Temp;
[m3,n3]=size(Aanvoertempvbuiten);
for i = 1:m3
    for j = 1:n3
        Aanvoertempvbuiten(i,j)=Gemiddeldewaardeperhalfuur(1,j)+Maxaanvoer+(Temp(i,j)*((Minaanvoer-Maxaanvoer)/(Maxbuiten-Minbuiten)));
        if Aanvoertempvbuiten(i,j)<55
            Aanvoertempvbuiten(i,j)=55;
        end
        if Aanvoertempvbuiten(i,j)>70
            Aanvoertempvbuiten(i,j)=70;
        end
    end
end
%COP ifv Aanvoertemperatuur
[m4,n4]=size(Aanvoertempvbuiten);
COPwaarden=zeros(m,n);
for i = 1:m4
    for j= 1:n4
        COPwaarden(i,j)=b(1)*Aanvoertempvbuiten(i,j)+b(2);
    end
end
%
dcechtverbruik=zeros(365,24);
for i = 1:365
    for j = 1:24
        dcechtverbruik(i,j)=Echtverbruik(i,2*j)+Echtverbruik(i,2*j-1);
    end
end

maxEV=max(max(dcechtverbruik));
minEV=min(min(dcechtverbruik));
Aantalvaluesingrafiekwarmtevraag=round((maxEV-minEV)*100,0)+1;
q=0; % variabele voor volgende for lus
Grafiekwarmtevraag=zeros(Aantalvaluesingrafiekwarmtevraag,2);
for i = 1:Aantalvaluesingrafiekwarmtevraag
    for j = 1:365
        for k = 1:24
            if dcechtverbruik(j,k)>= q && dcechtverbruik(j,k) < (minEV+i/100)
                Grafiekwarmtevraag(i,2)=Grafiekwarmtevraag(i,2)+1;
            end
        end
    end
    q=minEV+i/100;
    Grafiekwarmtevraag(i,1)=(minEV+(i-1)/100)*Aantalapp; %warmtevraag voor het totale gebouw
end

[m5,~]=size(Grafiekwarmtevraag);

for i=1:(m5-1)
    Grafiekwarmtevraag(m5-i,2)=Grafiekwarmtevraag(m5-i,2)+Grafiekwarmtevraag((m5-(i-1)),2);
end

%maken van elektriciteitsvraag
Elektriciteitsvraag=zeros(m,n);
for i = 1:m
    for j = 1:n
        Elektriciteitsvraag(i,j)=Echtverbruik(i,j)/COPwaarden(i,j);
    end
end
%duration curve elektriciteitsvraag
dcelekvraag=zeros(365,24);
for i = 1:365
    for j = 1:24
        dcelekvraag(i,j)=Elektriciteitsvraag(i,2*j)+Elektriciteitsvraag(i,2*j-1);
    end
end
Maxelek=max(max(dcelekvraag));
Minelek=min(min(dcelekvraag));
Aantalvaluesingrafiekelektriciteitsvraag=round((Maxelek-Minelek)*100,0)+1;
q=0; % variabele voor volgende for lus
Grafiekelektriciteitsvraag=zeros(Aantalvaluesingrafiekelektriciteitsvraag,2);
for i = 1:Aantalvaluesingrafiekelektriciteitsvraag
    for j = 1:365
        for k = 1:24
            if dcelekvraag(j,k)>= q && dcelekvraag(j,k) < (Minelek+i/100)
                Grafiekelektriciteitsvraag(i,2)=Grafiekelektriciteitsvraag(i,2)+1;
            end
        end
    end
    q=Minelek+i/100;
    Grafiekelektriciteitsvraag(i,1)=(Minelek+(i-1)/100)*Aantalapp; %warmtevraag voor het totale gebouw
end
[m5,n5]=size(Grafiekelektriciteitsvraag);
GE=Grafiekelektriciteitsvraag; % zelfde tabel voor economisch model

for i=1:(m5-1)
    Grafiekelektriciteitsvraag(m5-i,2)=Grafiekelektriciteitsvraag(m5-i,2)+Grafiekelektriciteitsvraag((m5-(i-1)),2);
end

% Economisch model zonder zonneboiler en buffervat
[m5,~]=size(Grafiekelektriciteitsvraag);
Overschotelektriciteit=0;
Overelektabel=zeros(m5,2);
for i=1:m5
    if Grafiekelektriciteitsvraag(i,1)<=warmtepompinout && waterradgemoutput>=warmtepompinout
        Overelektabel(i,1)=warmtepompinout-Grafiekelektriciteitsvraag(i,1);
        Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabel(i,1)=Overelektabel(i,1)+waterradgemoutput-warmtepompinout;
        Overschotelektriciteit=Overschotelektriciteit+Overelektabel(i,1)*GE(i,2);
    elseif Grafiekelektriciteitsvraag(i,1)>warmtepompinout && waterradgemoutput>=warmtepompinout
        Overelektabel(i,1)=Overelektabel(i,1)+waterradgemoutput-warmtepompinout;
        Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overschotelektriciteit=Overschotelektriciteit+Overelektabel(i,1)*GE(i,2);
    elseif waterradgemoutput<warmtepompinout
        if waterradgemoutput>=Grafiekelektriciteitsvraag(i,1)
            Overelektabel(i,1)=waterradgemoutput-Grafiekelektriciteitsvraag(i,1);
            Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2);
            Overschotelektriciteit=Overschotelektriciteit+Overelektabel(i,1)*GE(i,2);
        else
            Overelektabel(i,1)=0;
            Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2);
        end
        
    end
end

% PLUS/MIN afhankelijk van waterrad
Overelektabelmin=Overelektabel;
Overschotelektriciteitmin=0;
watermin=waterradgemoutput-waterradgemoutput*waterradplusmin;
for i=1:m5
    if Grafiekelektriciteitsvraag(i,1)<=warmtepompinout && watermin>=warmtepompinout
        Overelektabelmin(i,1)=warmtepompinout-Grafiekelektriciteitsvraag(i,1);
        Overelektabelmin(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabelmin(i,1)=Overelektabelmin(i,1)+watermin-warmtepompinout;
        Overschotelektriciteitmin=Overschotelektriciteitmin+Overelektabelmin(i,1)*GE(i,2);
    elseif Grafiekelektriciteitsvraag(i,1)>warmtepompinout && watermin>=warmtepompinout
        Overelektabelmin(i,1)=0;
        Overelektabelmin(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabelmin(i,1)=Overelektabelmin(i,1)+watermin-warmtepompinout;
        Overschotelektriciteitmin=Overschotelektriciteitmin+Overelektabelmin(i,1)*GE(i,2);
    elseif watermin<warmtepompinout
        if watermin>=Grafiekelektriciteitsvraag(i,1)
            Overelektabelmin(i,1)=watermin-Grafiekelektriciteitsvraag(i,1);
            Overschotelektriciteitmin=Overschotelektriciteitmin+Overelektabelmin(i,1)*GE(i,2);
        else
            Overelektabelmin(i,1)=0;
            Overelektabelmin(i,2)=Grafiekelektriciteitsvraag(i,2);
        end
    end
end
Overelektabelmax=Overelektabel;
Overschotelektriciteitmax=0;
watermax=waterradgemoutput+waterradgemoutput*waterradplusmin;
for i=1:m5
    if Grafiekelektriciteitsvraag(i,1)<=warmtepompinout && watermax>=warmtepompinout
        Overelektabelmax(i,1)=warmtepompinout-Grafiekelektriciteitsvraag(i,1);
        Overelektabelmax(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabelmax(i,1)=Overelektabelmax(i,1)+(watermax-warmtepompinout);
        Overschotelektriciteitmax=Overschotelektriciteitmax+Overelektabelmax(i,1)*GE(i,2);
    elseif Grafiekelektriciteitsvraag(i,1)>warmtepompinout && watermax>=warmtepompinout
        Overelektabelmax(i,1)=0;
        Overelektabelmax(i,1)=Overelektabelmax(i,1)+(watermax-warmtepompinout);
        Overschotelektriciteitmax=Overschotelektriciteitmax+Overelektabelmax(i,1)*GE(i,2);
        
    elseif watermax<warmtepompinout
        if watermax>=Grafiekelektriciteitsvraag(i,1)
            Overelektabelmax(i,1)=watermax-Grafiekelektriciteitsvraag(i,1);
            Overelektabelmax(i,2)=Grafiekelektriciteitsvraag(i,2);
        else
            Overelektabelmax(i,1)=0;
            Overelektabelmax(i,2)=Grafiekelektriciteitsvraag(i,2);
        end
    end
end
%Inkomsten
kwhovertomoney=[ watermin Overschotelektriciteitmin*exportgas;waterradgemoutput Overschotelektriciteit*exportgas;watermax Overschotelektriciteitmax*exportgas];
%Variabele kosten
extragas=0;
for i = 1:365
    for j = 1:24
        if dcelekvraag(i,j)*Aantalapp>warmtepompinout
            extragas=extragas+(dcelekvraag(i,j)*Aantalapp-warmtepompinout)*(COPwaarden(i,2*j)+COPwaarden(i,2*j-1))/2*Prijsgas;
        end
    end
end
variabelekostperjaar= Overschotelektriciteit*exportgas+Verwarmingskost-extragas;
Jaar=1:20;

%Piekduur
Piekelek=dcelekvraag*Aantalapp;
for i = 1:365
    for j = 1:24
        if dcelekvraag(i,j)*Aantalapp<=warmtepompinout
            Piekelek(i,j)=0;
        end
    end
end
%hoe lang duren die pieken
Piekduur=zeros(m,3);
kollominpiekduur=1;
for i = 1:365-1
    for j=1:24-1
        if Piekelek(i,j+1)>0 && Piekelek(i,j)>0
            Piekduur(i,kollominpiekduur)=Piekduur(i,kollominpiekduur)+0.5;
        elseif Piekelek(i,j+1)==0 && Piekelek(i,j)>0
            kollominpiekduur=kollominpiekduur+1;
        end
    end
    kollominpiekduur=1;
end
% duration curve van de pieken
maxpiekduur=max(max(Piekduur));
minpiekduur=min(min(Piekduur));
durationpiek=zeros(maxpiekduur-minpiekduur,2);
for deltapiek=1:(maxpiekduur-minpiekduur)
    for i= 1:m
        for j=1:3
            if Piekduur(i,j)== deltapiek+minpiekduur
                durationpiek(deltapiek,2)=durationpiek(deltapiek,2)+0.5;
            end
            durationpiek(deltapiek,1)=minpiekduur+deltapiek;
        end
    end
end
for i=1:(maxpiekduur-minpiekduur-1)
    durationpiek(maxpiekduur-minpiekduur-i,2)=durationpiek(maxpiekduur-minpiekduur-i,2)+durationpiek((maxpiekduur-minpiekduur-(i-1)),2);
end

%Omrekenen van Elektriciteit naar Warmtepompwarmte
Warmtepompwarmte = ones(365,24);
for i = 1:24
    for j = 1:365
        Warmtepompwarmte(j,i)=Warmtepompvermogen(j,i)*cteCOP;
    end
end
%Veranderen van Warmtevraag per halfuur naar een matrix met uurlijke waarden
Warmtevraaginuur = ones(365,48);
for i = 1:2:48
    for j = 1:365
        Warmtevraaginuur(j,i) = Warmtevraag(j,i)+Warmtevraag(j,i+1);
    end
end
for i = 2:25
    Warmtevraaginuur(:,i)=[];
end
Warmtepomptekort = ones(365,24);
for i = 1:24
    for j = 1:365
        Warmtepomptekort(j,i)=Warmtevraaginuur(j,i)-Warmtepompwarmte(j,i);
        if(Warmtepomptekort(j,i)<0)
            Warmtepomptekort(j,i)=0;
        end
    end
end
Warmtepompteveel = ones(365,24);
for i = 1:24
    for j = 1:365
        Warmtepompteveel(j,i) = Warmtepompwarmte(j,i)-Warmtevraaginuur(j,i);
        if(Warmtepompteveel(j,i)<0)
            Warmtepompteveel(j,i) = 0;
        end
    end
end
%Bufferberekeningen
Bufferopslag = zeros(365,24);
Warmtevraagmetbuffer = zeros(365,24);
Bufferoverschot = zeros(365,24);
for j = 1:365
    for i = 1:24
        if(i==1 && j == 1) %Opstart van Bufferopslag
            Bufferopslag(j,i) = Warmtepompteveel(j,i);
        elseif (i>1) %Uren na opstart
            Bufferopslag(j,i) = Warmtepompteveel(j,i)+Bufferoverschot(j,i-1);
        elseif (i == 1 && j>1) %Dag verder
            Bufferopslag(j,i) = Bufferopslag(j-1,24);
        end
        
        if(Bufferopslag(j,i)>W_tank)
            Bufferopslag(j,i)=W_tank;
        end
        
        if(Warmtevraaginuur(j,i)>Warmtepompwarmte(j,i))%Ontladen bij hogere warmtevraag dan levering
            Warmtevraagmetbuffer(j,i) = Warmtevraaginuur(j,i)-Bufferopslag(j,i); %Volledige bufferdepletie
            if(Warmtevraagmetbuffer(j,i)<Warmtepompwarmte(j,i)) %Bufferdepletie als niet volledige buffer leeg geraakt
                Bufferoverschot(j,i) = Warmtepompwarmte(j,i)-Warmtevraagmetbuffer(j,i); %Berekening van overblijvende buffer
                Warmtevraagmetbuffer(j,i) = Warmtepompwarmte(j,i); %Gecorrigeerde warmtevraag met buffer
            end
        elseif(Warmtevraaginuur(j,i)<=Warmtepompwarmte(j,i)) %Opladen Buffer
            Warmtevraagmetbuffer(j,i) = Warmtevraaginuur(j,i); %Ogenblikkelijke warmtevraag
            Bufferoverschot(j,i) = Bufferopslag(j,i);
        end
    end
end

%% Plot Data
plot([],[])
%% Name GUI Window
f.Name = 'Dijlenmolens GUI';
%% Position GUI Window
movegui(f,'center');
%% Make GUI visible
f.Visible = 'on';

%% Button Callback
    function hbutton1_Callback(source,eventdata)
        % Display plot of the currently selected data.
        cla reset
        plot(Aanvoertemp,COP);
        grid on
        hold on
        b=polyfit(Aanvoertemp,COP,1); %y=ax+b -> b(1) = a, b(2) = b
        X=1:80;
        plot(X,b(1)*X+b(2))
        xlabel('aanvoertemperatuur (°C)')
        ylabel('COP')
        title('Stooklijn')
    end
    function hbutton2_Callback(source,eventdata)
        % Display plot of the currently selected data.
        cla reset
        plot(Grafiekwarmtevraag(:,2),Grafiekwarmtevraag(:,1));
        grid on
        xlabel('Tijd (uur)')
        ylabel('Warmtevraag (kWh)')
        title('Warmtevraag')
    end
    function hbutton3_Callback(source,eventdata)
        cla reset
        plot(Grafiekelektriciteitsvraag(:,2),Grafiekelektriciteitsvraag(:,1))
        grid on
        xlabel('Tijd (uur)')
        ylabel('Elektriciteitsvraag (kWh)')
    end
    function hbutton4_Callback(source,eventdata)
        cla reset
        plot(Grafiekelektriciteitsvraag(:,2),Grafiekelektriciteitsvraag(:,1));
        grid on
        hold on
        plot(Overelektabel(:,2),Overelektabel(:,1))
        hold on
        plot(Overelektabelmin(:,2),Overelektabelmin(:,1))
        hold on
        plot(Overelektabelmax(:,2),Overelektabelmax(:,1))
        xlabel(' Tijd (uren)')
        ylabel(' Elektriciteitsvraag/overschot (kWh)')
        title('Elektriciteitsoverschot voor export afhankelijk van productie waterrad')
    end
    function hbutton5_Callback(source,eventdata)
        cla reset
        plot(kwhovertomoney(:,1),kwhovertomoney(:,2));
        axis([10 30 0 10000])
        grid on
        xlabel('gemiddelde productie waterrad (kWh)')
        ylabel('export inkomsten (euro)')
        title('export inkomsten')
    end
    function hbutton6_Callback(source,eventdate)
        cla reset
        plot(Jaar,variabelekostperjaar*Jaar);
        hold on
        plot(Jaar,ones(size(Jaar))*totalekost)
        grid on
        xlabel('jaren')
        ylabel('prijs')
        title('terugverdientijd')
    end
    function hbutton7_Callback(source,eventdate)
        cla reset
        plot(durationpiek(:,2),durationpiek(:,1));
        grid on
        xlabel('aantal keren per jaar dat warmtepomp onvoldoende is')
        ylabel('lengte van de piek (uur) ')
        title('piekduur')
        axis([0 90 0 10])
    end
    function hbutton8_Callback(source,eventdate)
        cla reset
        Dag = 1; %Starting Values
        x_time = time(Dag,:);
        y_1 = Warmtevraaginuur(Dag,:);
        y_2 = Warmtepompwarmte(Dag,:);
        y_3 = Warmtevraagmetbuffer(Dag,:);
        y_4 = Bufferopslag(Dag,:);
        %New UI Generation
        f2 = figure('Visible','on','units','normalized','outerposition',[0,0,1,1],'Numbertitle','off');
        hBufferbutton1 = uicontrol('Style','pushbutton',...
            'String','Dag Verder','Position',[900,750,90,25],...
            'Callback',@hBufferbutton1_Callback);
        hBufferbutton2 = uicontrol('Style','pushbutton',...
            'String','Dag Terug','Position',[650,750,90,25],...
            'Callback',@hBufferbutton2_Callback);
        hMenuPopup = uicontrol('Style', 'popup',...
            'String', {datestr(Datevector)},...
            'Position', [750 725 140 50],...
            'Callback', @hMenuPopup_Callback);
        %Plotting of Buffer Graphs
        ax1 = subplot(2,1,1);
        plot(x_time,y_1)
        title(datestr(Datevector(Dag)));
        hold on
        plot(x_time,y_2)
        plot(x_time,y_3)
        legend('Warmtevraag','Warmtepompgeleverde Warmte','Gebufferde Warmtevraag')
        xtickformat('hh:mm')
        xticks(0:hours(1):hours(23))
        xlabel('Tijd')
        ylabel('Verbruik (kWh)')
        ax2 = subplot(2,1,2); %Subplot van de Bufferopslag in kWh
        plot(x_time,y_4)
        xlabel('Tijd')
        ylabel('Bufferopslag (kWh)')
        xtickformat('hh:mm')
        xticks(0:hours(1):hours(23))
        function hBufferbutton1_Callback(source,eventdata)
            cla(ax1); cla(ax2);
            if(Dag == 365)
                Dag = 1;
            else
                Dag = Dag + 1;
            end
            set(hMenuPopup,'value',Dag);
            x_time = time(Dag,:);
            y_1 = Warmtevraaginuur(Dag,:);
            y_2 = Warmtepompwarmte(Dag,:);
            y_3 = Warmtevraagmetbuffer(Dag,:);
            y_4 = Bufferopslag(Dag,:);
            subplot(2,1,1)
            plot(x_time,y_1)
            hold on
            plot(x_time,y_2)
            plot(x_time,y_3)
            legend('Warmtevraag','Warmtepompgeleverde Warmte','Gebufferde Warmtevraag')
            xtickformat('hh:mm')
            xticks(0:hours(1):hours(23))
            xlabel('Tijd')
            ylabel('Verbruik (kWh)')
            subplot(2,1,2) %Subplot van de Bufferopslag in kWh
            plot(x_time,y_4)
            xlabel('Tijd')
            ylabel('Bufferopslag (kWh)')
            xtickformat('hh:mm')
            xticks(0:hours(1):hours(23))
        end
        function hBufferbutton2_Callback(source,eventdata)
            cla(ax1); cla(ax2);
            if(Dag == 1)
                Dag = 365;
            else
                Dag = Dag - 1;
            end
            set(hMenuPopup,'value',Dag);
            x_time = time(Dag,:);
            y_1 = Warmtevraaginuur(Dag,:);
            y_2 = Warmtepompwarmte(Dag,:);
            y_3 = Warmtevraagmetbuffer(Dag,:);
            y_4 = Bufferopslag(Dag,:);
            subplot(2,1,1)
            plot(x_time,y_1)
            hold on
            plot(x_time,y_2)
            plot(x_time,y_3)
            legend('Warmtevraag','Warmtepompgeleverde Warmte','Gebufferde Warmtevraag')
            xtickformat('hh:mm')
            xticks(0:hours(1):hours(23))
            xlabel('Tijd')
            ylabel('Verbruik (kWh)')
            subplot(2,1,2) %Subplot van de Bufferopslag in kWh
            plot(x_time,y_4)
            xlabel('Tijd')
            ylabel('Bufferopslag (kWh)')
            xtickformat('hh:mm')
            xticks(0:hours(1):hours(23))
        end
        function hMenuPopup_Callback(source,eventdata)
            cla(ax1); cla(ax2);
            selected_date = source.Value;
            Dag = selected_date;
            x_time = time(Dag,:);
            y_1 = Warmtevraaginuur(Dag,:);
            y_2 = Warmtepompwarmte(Dag,:);
            y_3 = Warmtevraagmetbuffer(Dag,:);
            y_4 = Bufferopslag(Dag,:);
            subplot(2,1,1)
            plot(x_time,y_1)
            hold on
            plot(x_time,y_2)
            plot(x_time,y_3)
            legend('Warmtevraag','Warmtepompgeleverde Warmte','Gebufferde Warmtevraag')
            xtickformat('hh:mm')
            xticks(0:hours(1):hours(23))
            xlabel('Tijd')
            ylabel('Verbruik (kWh)')
            subplot(2,1,2) %Subplot van de Bufferopslag in kWh
            plot(x_time,y_4)
            xlabel('Tijd')
            ylabel('Bufferopslag (kWh)')
            xtickformat('hh:mm')
            xticks(0:hours(1):hours(23))
        end
    end
end