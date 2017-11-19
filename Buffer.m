clc; clear;
c_p = 4.18; %kJ/kgK warmtecoefficient van water

fillrate = 0.5; %Liter/s
V_tank = 2000; %Tank Volume in Liter

t_filled = V_tank/fillrate; %tijd tot de tank volledig gevuld is
t=0:1:t_filled;
water_volume_in_tank = fillrate*t;

T_storage_winter = 273+70;
T_discharged_winter = 273+50;

Warmtevraaginhalfuur = xlsread('36AppWarmtevraag');
Warmtepompvermogen = 18*ones(365,24); %Productieprofiel van warmtepomp, deze gaat nog nauwkeuriger worden met waterradproductieprofiel
cteCOP = 3.6; %berekende sCOP
W_water = water_volume_in_tank*c_p*(T_storage_winter-T_discharged_winter)*(1/3600); %opgeslagen warmte in kWh
W_tank = max(W_water); %maximaal opgslaagbare warmte in de tank

firstDay = datetime(2016,1,1);
lastDay = datetime(2016,12,31);
Datevector = firstDay:caldays(1):lastDay;

time=ones(365,24);
for days = 1:365
    time(days,:) = 0:1:23;
end
time = hours(time);

Warmtepompwarmte = ones(365,24);
for i = 1:24
    for j = 1:365
           Warmtepompwarmte(j,i)=Warmtepompvermogen(j,i)*cteCOP;
    end
end

Warmtevraaginuur = turnHalfHourtoHourBySumming(Warmtevraaginhalfuur);

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
Dag = input('Geef dag van het jaar in nummer.  ');
x_time = time(Dag,:);
y_1 = Warmtevraaginuur(Dag,:);
y_2 = Warmtepompwarmte(Dag,:);
y_3 = Warmtevraagmetbuffer(Dag,:);
y_4 = Bufferopslag(Dag,:);
f = figure('Visible','on','units','normalized','outerposition',[0,0,1,1],'Numbertitle','off');
subplot(2,1,1)
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

subplot(2,1,2) %Subplot van de Bufferopslag in kWh
plot(x_time,y_4)
xlabel('Tijd')
ylabel('Bufferopslag (kWh)')
xtickformat('hh:mm')
xticks(0:hours(1):hours(23))