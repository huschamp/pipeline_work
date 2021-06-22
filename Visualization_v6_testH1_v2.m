clc
clear all
Table=readtable('maxpik838.xlsx');
%Table=readtable('minpik864.xlsx');
timein=table2array(Table(:,8));
timein=timein(1:(length(timein)-1));
timeout=table2array(Table(:,15));

datain=table2array(Table(:,3:6));
datain=datain(1:(length(datain)-1),:);
dataout=table2array(Table(:,10:13));
n=4 %число кранов
%n=3;
datain=round(datain*10000)/10000;
dataout=round(dataout*10000)/10000;
%
%datain=dataout;
%timein=timeout;
sized=size(datain);
sized=sized(2);
if sized>n
    for i=1:(sized-n)
        datain(:,sized)=[];
    end;
end;
%

CrivX=[];
CrivY=[];
MAX=[];
T_MAX=[];
for h=1:n
    %h=4;
    datainR=[datain(:,h) timein];
    dataoutR=[dataout(:,h) timeout];

    datain_desc=sortrows(datainR,1);
    dataout_desc=sortrows(dataoutR,1);
    for i=length(datain_desc(:,1)):-1:1
        if isnan(datain_desc(i,1))==1
            datain_desc(i,:)=[];
        end;
    end;
    datain_desc=datain_desc(end:-1:1,:);
    
    k=6;
    datainready=[];
    znin=sort(unique(datain_desc(:,1)),'descend');
    d1=sortrows(datain_desc,2);
    for i=1:k
        if i==1
            index=find(d1(:,1)==max(d1(:,1)));
            datainready=d1(index,:)%datain_desc(i,:);
        else
            if length(datainready(:,1))==1
                j=find(d1(:,1)==datainready(:,1));
                %j=length(j);
                if length(j)>1
                    k=k+1;
                    datainready=[datainready;d1(j(length(j)),:)];
                else             
                    jj=find(d1(:,1)==znin(i));
                    datainready=[datainready;d1(jj(1),:)];
                end;
            else
                jj=find(d1(:,1)==znin(i));
                datainready=[datainready;d1(jj(1),:)];
            end;
        end;
    end;
    
    datainready=sortrows(datainready,2);
    din=datainready;
    k=length(din(:,1)); %%%%%%%%%%%%%%%%%%%%%%??
    [maxz maxi]=max(datainready(:,1));
    
    index=find(din(:,1)==maxz);
    index(1)=[];
    if length(index)>0
        index(length(index))=[];
    end;
    for i=1:length(din(:,1))
        if length(find(index==i))>0
            din(i,:)=[];
        end;
    end;
    k=length(din(:,1)); %%%%%%%%%%%%%%%%%%%%%%??
    
    for i=maxi:k
        if i<=k-2
            y1=din(i,1);
            y2=din(i+2,1);
            x1=din(i,2);
            x2=din(i+2,2);
            x=din(i+1,2);
            y=y1+((y2-y1)/(x2-x1))*(x-x1);
            if din(i+1,1)<y
                din(k,:)=[];
                k=k-1;
            end;
        end;
    end;
    for i=maxi:-1:1  %maxi-1
        if i>=3
            y1=din(i,1);
            y2=din(i-2,1);
            x1=din(i,2);
            x2=din(i-2,2);
            x=din(i-1,2);
            y=y1+((y2-y1)/(x2-x1))*(x-x1);
            if din(i-1,1)<y
                k=k-1;
                din(1,:)=[];
            end;
        end;
    end;
    [maxz maxi]=max(din(:,1));
    %%%%%OSNOVNAY_CHAST
    if h==1
        din(1,:)=[];
        k=k-1;
        maxi=maxi-1;
    end;
    z1=0;
    z2=0;
%      x=din(:,2);
%      xx=linspace(min(x),max(x),1000);
%      y=din(:,1);
%      sp=interp1(x,y,xx,'nearest')
%      plot(xx,sp,'r-',x,y,'b.')
     if h==1
         eps=[0.015 0.001 0.015];
     end;
    if (h==1) && (n==3)
        eps=[0.025 0.002 0.02];
    end;
    if h==2
        eps=[0.012 0.0005 0.008];
    end;
    if h==3
        eps=[0.012 0.0005 0.008];
    end;
    if h==4
        eps=[0.009 0.0005 0.008];
    end;
    
        t=0:0.01:1;
        rangey=din(maxi,1):0.04:din(maxi,1)+15; %3 и 4 норм
        if h==1
            rangex=(din(1,2)-z1):12:(din(k,2)+165-z2);
        else
            rangex=(din(1,2)-z1):12:(din(k,2)+50-z2);
        end;
        if (h==1) && (length(timein)==length(timeout))
            rangex=(din(1,2)-z1):12:(din(k,2)-z2);
        end;
        Yin=[];
        Xin=[];

        for i=1:length(rangex)
            for j=1:length(rangey)
                FinY=(1-t).^2.*din(1,1)+2.*t.*(1-t).*rangey(j) +t.^2.*din(k,1);
                FinX=(1-t).^2.*(din(1,2)-z1)+2.*t.*(1-t).*rangex(i) +t.^2.*(din(k,2)-z2);
                Yin=[Yin;FinY];
                Xin=[Xin;FinX];
            end;
        end;
        L=length(Xin(:,1));

        dev=10;
        s=[];
        M=[];
        for g=1:L
            if (max(Yin(g,:))>max(din(:,1))+eps(1)) || (max(Yin(g,:))<max(din(:,1))-eps(1))
                dev=10;
                s=[s dev];
                continue;
            end;
            
%             [maxg maxgindex]=max(Yin(g,:));
%             maxdin=max(din(:,1));
%             maxdinindex=find(din(:,1)==maxdin);
%             mleft=find(Xin(g,:)<=din(maxdinindex(1)-1,2));
%             mleft=mleft(length(mleft));
%             mright=find(Xin(g,:)<=din(maxdinindex(end)+1,2));
%             mright=mright(length(mright));
%             if (maxgindex<=mleft) || (maxgindex>=mright) 
%                 dev=10;
%                 s=[s dev];
%                 continue;
%             end;
            %f=find(d1(:,1)==din(i-1,1));
            [maxg maxgindex]=max(Yin(g,:));
            maxdin=max(din(:,1));
            maxdinindex=find(din(:,1)==maxdin);
            f=find(d1(:,1)==din(maxdinindex(1)-1,1));
            mleft=find(Xin(g,:)<=d1(f(length(f)),2));
            mleft=mleft(length(mleft));
            mright=find(Xin(g,:)<=din(maxdinindex(end)+1,2));
            mright=mright(length(mright));
            if (maxgindex<=mleft) || (maxgindex>=mright) 
                dev=10;
                s=[s dev];
                continue;
            end;
            
            err=0;
            M=[];
            for i=2:k
                j2=find(Xin(g,:)<=din(i,2));
                j2=j2(length(j2));

                f=find(d1(:,1)==din(i-1,1));
                i1=find(Xin(g,:)>=d1(f(1),2));
                i1=i1(1);
                j1=find(Xin(g,:)<=d1(f(length(f)),2));
                j1=j1(length(j1));
                
                if din(i,1)>din(i-1,1)
                    if length(find(Yin(g,i1:j1)>(din(i,1)+eps(2))))>0
                        dev=10;
                        err=1;
                        s=[s dev];
                        break;
                    end;
                end;
                epsilon=0.015;
                if din(i-1,1)~=din(i,1)
                    if (length(find(Yin(g,j1:j2)>(din(i,1)+eps(3))))==length(Yin(g,j1:j2))) || (length(find(Yin(g,j1:j2)<(din(i,1)-eps(3))))==length(Yin(g,j1:j2)))
                        dev=10;
                        err=1;
                        s=[s dev];
                        break;
                    end;
                end;
                M=[M abs(min(Yin(g,j2)-din(i,1)))];
            end;
            if err==0
                dev=max(M);
                s=[s dev];
            end;
        end;
    [OPT i]=min(s);
    interest=abs(max(Xin(i,:))-max(din(:,2)));
    R=[];
    if h==1
        fin=find(s~=10);
        for j=1:length(find(s~=10))
            R=[R abs(max(Xin(fin(j),:))-max(din(:,2)))];
        end;
        [OPTimalX i]=min(R);
        i=fin(5);
    end;
    MAX=[MAX max(Yin(i,:))];
    find(Yin(i,:)==max(Yin(i,:)))
    T_MAX=[T_MAX Xin(i,find(Yin(i,:)==max(Yin(i,:))))];
    CrivX=[CrivX;Xin(i,:)];
    CrivY=[CrivY;Yin(i,:)];
    plot(d1(:,2),d1(:,1),'.',Xin(i,:),Yin(i,:))
    figure
end;
plot(timein,datain,'.', CrivX',CrivY')