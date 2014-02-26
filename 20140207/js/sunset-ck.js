/*
 JavaScript Sun Table Calculator
© 2001 Juergen Giesen
http://www.jgiesen.de
*/function setMyMenu(e,t,n,r,i){}function currentTime(){var e=new Date,t=e.getTimezoneOffset();t>=1380&&(t-=1440);var n=-t/60;utHours=e.getHours()-n;LatLong=document.myform.location.options[document.myform.location.selectedIndex].value;document.myform.location.selectedIndex!=0&&getLocation(LatLong);var r=2*locOffset;locOffset>0&&(r-=1);document.myform.timezone.options[Math.abs(r)].selected=!0;hours=utHours+locOffset;hours<0&&(hours+=24);hours>23&&(hours-=24);minutes=e.getMinutes();seconds=e.getSeconds();seconds=0;year=e.getYear();year<1900&&(year+=1900);month=e.getMonth()+1;day=e.getDate();utMinutes=minutes;utMonth=month;utYear=year;document.myform.localDateText.value=writeDateTime(hours,day,month);setMyMenu(day,month,year,hours,minutes)}function writeDateTime(e,t,n){var r=new Date(year,n-1,t,hours,minutes,seconds),i=dayName[r.getDay()],s=" "+year+", "+i+" "+document.myform.Monat.options[n-1].value+" ";t<10?s=s+"0"+t:s+=t;s+=" at ";e<10?s=s+"0"+e+":":s=s+e+":";minutes<10?s=s+"0"+minutes:s+=minutes;return s}function getDateTime(){var e=getLocation("52.38/4.89/+1*Amsterdam"),t=e.length,n="",r=0,i=0,s=0;for(var o=0;o<t;o++)if(e.charAt(o)=="/"){r=o;break}for(var o=r+1;o<t;o++)if(e.charAt(o)=="/"){i=o;break}for(var o=i+1;o<t;o++)if(e.charAt(o)=="*"){s=o;break}locName=e.substring(s+1,t);offset=-60*document.myform.timezone.options[document.myform.timezone.selectedIndex].value;locOffset=-offset/60;var u=new Date;hours=u.getHours();minutes=u.getMinutes();seconds=u.getSeconds();year=u.getYear();year<1900&&(year+=1900);month=u.getMonth()+1;day=u.getDate();var a=u.getTime();a-=locOffset*3600*1e3;var f=new Date(a);utHours=f.getHours();utDay=f.getDate();utMonth=f.getMonth()+1;utYear=f.getYear();utYear<1900&&(utYear+=1900);seconds=0;UT=utHours+minutes/60+seconds/3600;document.myform.localDateText.value=writeDateTime(hours,day,month);document.myform.UTimeText.value=writeDateTime(utHours,utDay,utMonth);setMyMenu(day,month,year,hours,minutes);EOT=eot(utDay,utMonth,utYear,UT);theJulDay()}function getTZone(){var e,t,n,r,i=document.myform.timezone.selectedIndex,s=document.myform.timezone.options[i].value,o=new Date(Date.UTC(year,month-1,utDay,utHours,minutes,seconds)),u=Date.UTC(year,month-1,day,utHours,minutes,seconds),a=o.getTimezoneOffset();a>=1380&&(a-=1440);var f=-a/60;u-=(f-s)*3600*1e3;var l=new Date(u);hours=l.getHours();minutes=l.getMinutes();seconds=l.getSeconds();year=l.getYear();year<1900&&(year+=1900);month=l.getMonth()+1;day=l.getDate();if(s==0){day=utDay;month=utMonth}if(s>0){if(month>utMonth&&Number(utHours)+Number(s)>=24){month=utMonth;day=utDay+1;daysInMonth(month,year);if(day>dIM){day=1;month+=1}}if(month==utMonth||Number(utHours)+Number(s)<24){day>utDay+1&&(day-=1);month>utMonth&&(day=1);if(Number(utHours)+Number(s)<24){day=utDay;month=utMonth;year=utYear}year>utYear&&day>1&&(day-=1)}}if(s<0){daysInMonth(month,year);if(day==dIM&&utDay==2){day=1;month=utMonth}if(day!=dIM&&utDay!=2){day<utDay-1&&(day+=1);if(month<utMonth){daysInMonth(month,year);day=dIM}if(Number(utHours)+Number(s)>=0){day=utDay;month=utMonth}year<utYear&&day<31&&(day+=1)}}setMyMenu(day,month,year,hours,minutes);document.myform.localDateText.value=writeDateTime(hours,day,month);calcRiseSet()}function ut(){utHours=hours-locOffset}function setUT(){seconds=0;var e=new Date(year,month-1,day,hours,minutes,seconds),t=e.getTime();t-=locOffset*3600*1e3;var n=new Date(t);utHours=n.getHours();utMinutes=n.getMinutes();utDay=n.getDate();utMonth=n.getMonth()+1;utYear=n.getYear();utYear<1900&&(utYear+=1900);UT=utHours+utMinutes/60+seconds/3600;document.myform.localDateText.value=writeDateTime(hours,day,month);document.myform.UTimeText.value=writeDateTime(utHours,utDay,utMonth);EOT=eot(utDay,utMonth,utYear,UT)}function getLocation(e){var t=e.length,n="",r=0,i=0,s=0;for(var o=0;o<t;o++)if(e.charAt(o)=="/"){r=o;break}for(var o=r+1;o<t;o++)if(e.charAt(o)=="/"){i=o;break}for(var o=i+1;o<t;o++)if(e.charAt(o)=="*"){s=o;break}var n=e.substring(0,r);lat=Number(n);var u=lat,a=e.substring(r+1,i);longit=Number(a);var f=longit,l=e.substring(i+1,s);locOffset=Number(l);var c=2*locOffset;locOffset>0&&(c-=1);locName=e.substring(s+1,t);var h=Date.UTC(year,utMonth-1,utDay,utHours,minutes,seconds);h+=locOffset*3600*1e3;var p=new Date,d=p.getTimezoneOffset();d>=1380&&(d-=1440);var v=-d/60;h-=v*3600*1e3;p=new Date(h);hours=p.getHours();day=p.getDate();month=p.getMonth()+1;year=p.getYear();year<1900&&(year+=1900);setMyMenu(day,month,year,hours,minutes);locOffset>=0?timezoneString="GMT + "+locOffset:timezoneString="GMT  "+locOffset;UT=utHours+minutes/60+seconds/3600;calculate();calcRiseSet()}function riseset(e,t,n,r){var i=Math.PI/180,s=Math.sin(i*twilight),o=declination(e,t,n,r),u=computeGHA(e,t,n,r);Y0=Math.sin(i*computeHeight(o,lat,longit,u))-s;o=declination(e,t,n,r+1);u=computeGHA(e,t,n,r+1);yPlus=Math.sin(i*computeHeight(o,lat,longit,u))-s;o=declination(e,t,n,r-1);u=computeGHA(e,t,n,r-1);yMinus=Math.sin(i*computeHeight(o,lat,longit,u))-s;ABOVE=yMinus>0;QUAD(yMinus,yPlus);switch(NZ){case 0:break;case 1:if(yMinus<0){UTRISE=r+zero1;RISE=!0}else{UTSET=r+zero1;SETT=!0}break;case 2:if(YE<0){UTRISE=r+zero2;UTSET=r+zero1}else{UTRISE=r+zero1;UTSET=r+zero2}RISE=!0;SETT=!0}}function QUAD(e,t){NZ=0;var n,r=.5*(e+t)-Y0,i=.5*(t-e),s=Y0;n=-i/(2*r);YE=(r*n+i)*n+s;var o=i*i-4*r*s;if(o>=0){DX=.5*Math.sqrt(o)/Math.abs(r);zero1=n-DX;zero2=n+DX;Math.abs(zero1)<=1&&(NZ+=1);Math.abs(zero2)<=1&&(NZ+=1);zero1<-1&&(zero1=zero2)}}function JulDay(e,t,n,r){n<1900&&(n+=1900);if(t<=2){t+=12;n-=1}A=Math.floor(n/100);JD=Math.floor(365.25*(n+4716))+Math.floor(30.6001*(t+1))+e-13-1524.5+r/24;return JD}function computeHeight(e,t,n,r){var i=Math.PI/180,s=t*i,o=e*i,u=Number(r)+Number(n);sinHeight=Math.sin(o)*Math.sin(s)+Math.cos(o)*Math.cos(s)*Math.cos(i*u);return Math.asin(sinHeight)/i}function declination(e,t,n,r){var i=Math.PI/180,s=JulDay(e,t,n,r),o=(s-2451545)/36525,u=280.46645+(36000.76983+3032e-7*o)*o,a=357.5291+(35999.0503-(1559e-7*o+4.8e-7*o)*o)*o;a=i*a;var f=(1.9146-.004817*o-14e-6*o*o)*Math.sin(a)+(.019993-101e-6*o)*Math.sin(2*a)+29e-5*Math.sin(3*a),l=u+f,c=125.04-1934.136*o,h=l-.00569-.00478*Math.sin(i*c),p=23.43929111111111-(46.815*o+59e-5*o*o-.001813*o*o*o)/3600,d=p+.00256*Math.cos(i*c),v=Math.sin(i*d)*Math.sin(i*h);v=Math.asin(v)/i;RA=Math.atan2(Math.cos(i*d)*Math.sin(i*h),Math.cos(i*h))/i;RA<0&&(RA+=360);return v}function computeGHA(e,t,n,r){var i=Math.PI/180,s,o,u,a;s=365*n+e+31*t-46;t<3?s+=Math.floor((n-1)/4):s=s-Math.floor(.4*t+2.3)+Math.floor(n/4);a=r/24;o=(a+s-722449)*.98564734+279.306;o*=i;u=-104.55*Math.sin(o)-429.266*Math.cos(o)+595.63*Math.sin(2*o)-2.283*Math.cos(2*o);u=u+4.6*Math.sin(3*o)+18.7333*Math.cos(3*o);u=u-13.2*Math.sin(4*o)-Math.cos(5*o)-Math.sin(5*o)/3+.5*Math.sin(6*o)+.231;u=u/240+360*(a+.5);u>360&&(u-=360);return u}function computeAzimut(e,t,n,r,i){var s=Math.PI/180,o,u,a,f=t*s,l=i*s;a=Math.cos(l)*Math.cos(f);o=(Math.sin(e*s)-Math.sin(f)*Math.sin(l))/a;u=Math.PI/2-Math.asin(o);u/=s;Math.sin(s*(Number(r)+Number(n)))<=0?u=u:u=360-u;return u}function GM_Sidereal_Time(e){var t,n,r,i;i=e-2400000.5;r=Math.floor(i);n=(i-r)*24;t=(r-51544.5)/36525;return 6.697374558+1.0027379093*n+(8640184.812866+(.093104-62e-7*t)*t)*t/3600}function LM_Sidereal_Time(e,t){var n=GM_Sidereal_Time(e),r=24*frac((n+t/15)/24);return HoursMinutesSeconds(r)}function sunL(e){var t=280.46645+36000.76983*e+3032e-7*e*e;t%=360;t<0&&(t+=360);return t}function deltaPSI(e){var t=Math.PI/180,n,r,i,s;i=sunL(e);s=218.3165+481267.8813*e;s%=360;s<0&&(s+=360);r=125.04452-1934.136261*e+.0020708*e*e+e*e*e/45e4;n=-17.2*Math.sin(t*r)-1.32*Math.sin(t*2*i)-.23*Math.sin(t*2*s)+.21*Math.sin(t*2*r);n/=3600;return n}function EPS(e){var t=Math.PI/180,n=sunL(e),r=218.3165+481267.8813*e,i=23.43929111111111-(46.815*e+59e-5*e*e-.001813*e*e*e)/3600,s=125.04452-1934.136261*e+.0020708*e*e+e*e*e/45e4,o=(9.2*Math.cos(t*s)+.57*Math.cos(t*2*n)+.1*Math.cos(t*2*r)-.09*Math.cos(t*2*s))/3600;return i+o}function eot(e,t,n,r){var i=Math.PI/180,s=(JulDay(e,t,n,r)-2451545)/36525,o=EPS(s),u=declination(e,t,n,r),a=sunL(s),f=deltaPSI(s),l=a-.0057183-RA+f*Math.cos(i*o);l>5&&(l-=360);l*=4;l=Math.round(1e3*l)/1e3;return l}function transit(e,t,n){var r=declination(e,t,n,utHours),i=12-eot(e,t,n,12)/60-longit/15,s=computeGHA(e,t,n,i),o=declination(e,t,n,i),u=computeHeight(o,lat,longit,s);maxElevString=Math.round(10*u)/10;Math.round(10*u)/10==Math.round(u)&&(maxElevString+=".0");i+=locOffset;transString=HoursMinutes(i)}function dayTable(e){var t,n,r,i,s,o,u,a,f,l,c,h,p,d,v,m;calcRiseSet();lat>=0?ns=" N":ns=" S";longit>=0?ew=" E":ew=" W";locOffset>=0?timezoneString="GMT + "+locOffset:timezoneString="GMT  "+locOffset;c=locName+",  "+timezoneString+"\n";document.myform.area.value+=c;i=c;c="Latitude: "+Math.abs(lat)+deg+ns+",  Longitude: "+Math.abs(longit)+deg+ew+"\n";document.myform.area.value+=c;s=c;document.myform.area.value+="Local time: "+document.myform.localDateText.value+"\n";o="Local time: "+document.myform.localDateText.value+"\n";c="Sunrise: "+document.myform.riseTime.value;c=c+"     Sunset: "+document.myform.setTime.value;c=c+"    Length of Day: "+lengthDayString+" h";document.myform.area.value+=c+"\n";u=c;c="Equation of Time: "+EOT+" min = "+HoursMinutes(EOT)+" min";document.myform.area.value+=c+"\n";a=c;var h=declination(utDay,utMonth,utYear,UT);h=Math.round(1e3*h)/1e3;c="Declination: "+h+" deg = "+DegreesMinutes(h);h>=0?c+=" N":c+=" S";document.myform.area.value+=c+"\n";f=c;transit(day,month,year);document.myform.area.value+="Transit at: "+transString+"   Altitude: "+maxElevString+deg+"\n";l="Transit at: "+transString+"&nbsp;&nbsp;&nbsp;Altitude: "+maxElevString+deg;document.myform.area.value+="Table of Altitude, Azimuth, and GHA Angles:\n";for(var g=-locOffset;g<=-locOffset+24;g++){m=g+locOffset;h=declination(utDay,utMonth,utYear,g);p=computeGHA(utDay,utMonth,utYear,g);d=computeHeight(h,lat,longit,p);v=computeAzimut(h,lat,longit,p,d);v=Math.round(10*v)/10;v<100&&(v="  "+v);v<10&&(v="   "+v);Math.round(v)==Math.round(10*v)/10&&(v+=".0");d=Math.round(10*d)/10;Math.round(d)==Math.round(10*d)/10&&(d+=".0");d<10&&d>0&&(d=" "+d);d>0&&(d="  "+d);d>-10&&d<0&&(d=" "+d);m<10?t=" "+m:t=m;t+=":00";n=d;r=v;str4=Math.round(10*p)/10;Math.round(p)==Math.round(10*p)/10&&(str4+=".0");if(e==1)document.myform.area.value+=t+"   "+n+"   "+r+"   "+str4+"\n";else{strA=new Array(t,n+deg,r+deg);strArray[m]=strA;text[0]=i;text[1]=s;text[2]=o;text[3]=u;text[4]=a;text[5]=f;text[6]=l}}e==2&&writeDayPage(strArray)}function table(e){var t="",n="",r;lat>=0?ns=" N":ns=" S";longit>=0?ew=" E":ew=" W";locOffset>=0?timezoneString="GMT + "+locOffset:timezoneString="GMT  "+locOffset;tDat=new Date;daysInMonth(month,year);var i="Sunrise, Sunset, Length of Day, Transit, and Max. Altitude\n";i=i+"Latitude: "+Math.abs(lat)+deg+ns+",  Longitude: "+Math.abs(longit)+deg+ew+",   "+locName+",  "+timezoneString+"\n";text[0]="Latitude: "+Math.abs(lat)+deg+ns+",  Longitude: "+Math.abs(longit)+deg+ew+",&nbsp;&nbsp;&nbsp;"+locName+",  "+timezoneString;document.myform.area.value+=i;var o=this.location.toString();if(o.lastIndexOf(String.fromCharCode(106,103,105,101,115))==-1&&o.lastIndexOf(string1)==-1){document.myform.area.value+=ok+"\n";dIM=1}for(var u=1;u<=dIM;u++){JD=JulDay(u,month,year,12);EOT=eot(u,month,year,12);day=u;RISE=!1;SETT=!1;twilight=-0.8333;for(var f=-locOffset;f<-locOffset+24;f++){riseset(u,month,year,f);if(RISE&&SETT)break}if(RISE||SETT){if(RISE){UTRISE+=locOffset;UTRISE>=24&&(UTRISE-=24);UTRISE<0&&(UTRISE+=24)}if(SETT){UTSET+=locOffset;UTSET>=24&&(UTSET-=24);UTSET<0&&(UTSET+=24)}lengthOfDay=UTSET-UTRISE}r=UTRISE-locOffset+hours-utHours;r<0&&(r+=24);r>=24&&(r-=24);RISE?str1=HoursMinutes(r):ABOVE?str1="visible":str1="--.--";r=UTSET-locOffset+hours-utHours;r<0&&(r+=24);r>=24&&(r-=24);if(SETT)str2=HoursMinutes(r);else if(ABOVE){str2="visible";lengthOfDay=24}else{str2="--.--";lengthOfDay=0}RISE=!1;SETT=!1;twilight=-6;for(var f=-locOffset;f<-locOffset+24;f++){riseset(u,month,year,f);if(RISE&&SETT)break}if(RISE||SETT){if(RISE){UTRISE+=locOffset;UTRISE>=24&&(UTRISE-=24);UTRISE<0&&(UTRISE+=24)}if(SETT){UTSET+=locOffset;UTSET>=24&&(UTSET-=24);UTSET<0&&(UTSET+=24)}}r=UTRISE-locOffset+hours-utHours;r<0&&(r+=24);r>=24&&(r-=24);RISE?t=HoursMinutes(r):ABOVE?t="visible":t="--.--";r=UTSET-locOffset+hours-utHours;r<0&&(r+=24);r>=24&&(r-=24);SETT?n=HoursMinutes(r):ABOVE?n="visible":n="--.--";st1=""+Math.round(100*lengthOfDay)/100;st2=""+Math.round(10*lengthOfDay)/10;a=Math.round(100*lengthOfDay)/100;Math.round(lengthOfDay)==a?i=a+".00 h":st1.length==st2.length?i=a+"0 h":i=a+" h";u<10?s="0":s="";transit(day,month,year);str3=transString;str4=maxElevString+deg;var l=declination(day,month,year,utHours);EOT=eot(day,month,year,12+locOffset);str5=DegreesMinutes(EOT);e==1&&(document.myform.area.value+=year+" "+monthName[month-1]+" "+s+day+"    "+str1+"    "+str2+"   "+i+"   "+str3+"  "+str4+"\n");strB=new Array(year,monthName[month-1],day,t,str1,str2,n,i,str3,str4,str5);strArray[u]=strB}e==2&&writeMonthPage(strArray,dIM)}function calcRiseSet(){offset=-60;locOffset=-offset/60;twilight=-0.8333;RISE=!1;SETT=!1;ok=OK.charAt(1)+OK.charAt(3)+OK.charAt(2)+OK.charAt(0);for(var e=-locOffset;e<-locOffset+24;e++){riseset(day,month,year,e);if(RISE&&SETT)break}if(RISE||SETT){if(RISE){hRise1=UTRISE+locOffset;hRise1>24&&(hRise1-=24);hRise1<0&&(hRise1+=24);var t=declination(utDay,utMonth,utYear,UTRISE),n=computeGHA(utDay,utMonth,utYear,UTRISE),r=computeHeight(t,lat,longit,n),i=computeAzimut(t,lat,longit,n,r);i=Math.round(10*i)/10;document.myform.azimRise.value=" "+i}if(SETT){hSet1=UTSET+locOffset;hSet1>24&&(hSet1-=24);hSet1<0&&(hSet1+=24);var t=declination(utDay,utMonth,utYear,UTSET),n=computeGHA(utDay,utMonth,utYear,UTSET),r=computeHeight(t,lat,longit,n),i=computeAzimut(t,lat,longit,n,r);i=Math.round(10*i)/10}lengthOfDay=hSet1-hRise1}var s=this.location.toString();if(s.lastIndexOf(String.fromCharCode(106,103,105,101,115))!=-1||s.lastIndexOf(string1)!=-1){RISE?str=HoursMinutes(hRise1):ABOVE?str="visible":str="--.--";document.myform.riseTime.value=" "+str;if(SETT)str=HoursMinutes(hSet1);else if(ABOVE){str="visible";lengthOfDay=24}else{str="--.--";lengthOfDay=0}}lengthOfDay=Math.round(100*lengthOfDay)/100;lengthDayString=lengthOfDay;twilight=document.myform.twilight.options[document.myform.twilight.selectedIndex].value;RISE=!1;SETT=!1;for(var e=-locOffset;e<-locOffset+24;e++){riseset(day,month,year,e);if(RISE&&SETT)break}if(RISE||SETT){if(RISE){hRise2=UTRISE+locOffset;hRise2>24&&(hRise2-=24);hRise2<0&&(hRise2+=24)}if(SETT){hSet2=UTSET+locOffset;hSet2>24&&(hSet2-=24);hSet2<0&&(hSet2+=24)}lengthOfDay=hSet2-hRise2}if(s.lastIndexOf(String.fromCharCode(106,103,105,101,115))!=-1||s.lastIndexOf(string1)!=-1){RISE?str=HoursMinutes(hRise2):ABOVE?str="visible":str="--.--";document.myform.tRiseTime.value=" "+str;if(SETT)str=HoursMinutes(hSet2);else if(ABOVE){str="visible";lengthOfDay=24}else{str="--.--";lengthOfDay=0}document.myform.tSetTime.value=" "+str}transit(day,month,year);document.myform.transText.value=" "+transString;document.myform.maxAltText.value=" "+maxElevString}function calculate(){var e,t,n=" "+Math.round(1e3*declination(utDay,utMonth,utYear,UT))/1e3;OK="ODMEDO";lat=e;longit=t;var r=declination(utDay,utMonth,utYear,UT),i=computeGHA(utDay,utMonth,utYear,UT),s=computeHeight(r,lat,longit,i);s=Math.round(10*s)/10;EOT=eot(utDay,utMonth,utYear,UT);n=" "+EOT;var o=computeAzimut(r,lat,longit,i,s);i=Math.round(10*i)/10;var u=RA/15;u>=24&&(u-=24)}function theDay(){daysInMonth(month,year);var e=document.myform.Tag.selectedIndex;if(e<dIM)day=e+1;else{day=dIM;document.myform.Tag.options[dIM-1].selected=!0}setUT()}function theMonth(){month=Number(document.myform.Monat.selectedIndex)+Number(1);theDay()}function theYear(){year=document.myform.Jahr.selectedIndex+1992;theDay()}function theHour(){hours=document.myform.Stunde.selectedIndex}function theMinute(){minutes=document.myform.Minute.selectedIndex;seconds=0}function theJulDay(){jd=JulDay(utDay,utMonth,utYear,UT);document.myform.JulDayText.value=" "+Math.round(1e4*jd)/1e4}function theDateTime(){theDay();theMonth();theYear();theHour();theMinute();setUT();theJulDay();calculate();calcRiseSet()}function getLatitude(){locName="User Input";str=document.myform.latitude.value;var e=str.length;for(var t=0;t<e;t++){c=str.charAt(t);if(c!="0"&&c!="1"&&c!="2"&&c!="3"&&c!="4"&&c!="5"&&c!="6"&&c!="7"&&c!="8"&&c!="9"&&c!="+"&&c!="."){console.log("Error on latitude value !\nEnter positive decimal degree value, e.g.:\n 52.34, and select North or South.");document.myform.latitude.value=0;break}}if(Math.abs(Number(str))>90){console.log("Latitude must less or equal to 90 degrees !");document.myform.latitude.value=0}lat=Number(str);document.myform.NorthSouth.selectedIndex==0?lat=Math.abs(lat):lat=-Math.abs(lat);lat>=0?ns=" N":ns=" S";document.myform.location.options[0].selected=!0}function getLongitude(){str=document.myform.longitude.value;var e=str.length;for(var t=0;t<e;t++){c=str.charAt(t);if(c!="0"&&c!="1"&&c!="2"&&c!="3"&&c!="4"&&c!="5"&&c!="6"&&c!="7"&&c!="8"&&c!="9"&&c!="+"&&c!="."){console.log("Error on longitude value !\nEnter positive decimal degree value, e.g.:\n 8.34 and select East or West.");document.myform.longitude.value=0;break}}if(Math.abs(Number(str))>180){console.log("Longitude must less or equal to 180 degrees !");document.myform.longitude.value=0}longit=Number(str);document.myform.EastWest.selectedIndex==0?longit=Math.abs(longit):longit=-Math.abs(longit);longit>=0?ew=" E":ew=" W";document.myform.location.options[0].selected=!0;calculate();calcRiseSet()}function getTwilight(){twilight=document.myform.twilight.options[document.myform.twilight.selectedIndex].value;theDateTime()}function clearTable(){document.myform.area.value=" "}function NS(){document.myform.NorthSouth.selectedIndex==0?lat=Math.abs(lat):lat=-Math.abs(lat);locName="User Input";document.myform.location.options[0].selected=!0;calcRiseSet()}function EW(){document.myform.EastWest.selectedIndex==0?longit=Math.abs(longit):longit=-Math.abs(longit);locName="User Input";document.myform.location.options[0].selected=!0;locName="User Input";calculate();calcRiseSet()}function yearTable(){var e=new Array(13),t=new Array(13),n;lat>=0?ns=" N":ns=" S";longit>=0?ew=" E":ew=" W";locOffset>=0?timezoneString="GMT + "+locOffset:timezoneString="GMT  "+locOffset;tDat=new Date;text[0]=year+",&nbsp;&nbsp;&nbsp;"+"Latitude: "+Math.abs(lat)+deg+ns+",  Longitude: "+Math.abs(longit)+deg+ew+",&nbsp;&nbsp;&nbsp;"+locName+",  "+timezoneString;twilight=-0.8333;for(var r=1;r<=31;r++){for(var i=1;i<=12;i++){daysInMonth(i,year);RISE=!1;SETT=!1;for(var s=-locOffset;s<-locOffset+24;s++){riseset(r,i,year,s);if(RISE&&SETT)break}if(RISE||SETT){if(RISE){UTRISE+=locOffset;UTRISE>=24&&(UTRISE-=24);UTRISE<0&&(UTRISE+=24)}if(SETT){UTSET+=locOffset;UTSET>=24&&(UTSET-=24);UTSET<0&&(UTSET+=24)}}n=UTRISE-locOffset+hours-utHours;n<0&&(n+=24);n>=24&&(n-=24);RISE?str1=HoursMinutes(n):ABOVE?str1="vis.":str1="--.--";n=UTSET-locOffset+hours-utHours;n<0&&(n+=24);n>=24&&(n-=24);SETT?str2=HoursMinutes(n):ABOVE?str2="vis.":str2="--.--";if(r<=dIM){e[i]=str1;t[i]=str2}else{e[i]="";t[i]=""}var o=this.location.toString();o.lastIndexOf(String.fromCharCode(106,103,105,101,115))==-1&&o.lastIndexOf(string1)==-1&&(t[i]=ok)}strB=new Array(r,e[1],t[1],e[2],t[2],e[3],t[3],e[4],t[4],e[5],t[5],e[6],t[6],e[7],t[7],e[8],t[8],e[9],t[9],e[10],t[10],e[11],t[11],e[12],t[12]);strArray[r]=strB}writeYearPage(strArray)}function HoursMinutes(e){var t=e;e=Math.abs(e);var n=Math.round(60*(e-Math.floor(e))),r;r=Math.floor(e);r<10&&(r="0"+r);n>=10?r=r+":"+n:r=r+":0"+n;n==60&&(r=Math.floor(e+1)+":00");return t<0?"-"+r:r}function HoursMinutesSeconds(e){var t=Math.floor(e),n=Math.floor(60*frac(e)),r=Math.round(60*(60*frac(e)-n)),i;n>=10?i=t+":"+n:i=t+":0"+n;r<10?i=i+":0"+r:i=i+":"+r;return" "+i}function DegreesMinutes(e){var t=e;e=Math.abs(e);var n=Math.round(60*(e-Math.floor(e))),r;n>=10?r=Math.floor(e)+deg+" "+n+"'":r=Math.floor(e)+deg+" 0"+n+"'";n==60&&(r=Math.floor(e+1)+":00"+"'");return t<0?"-"+r:r}function daysInMonth(e,t){var n=31;e-=1;if(e==0||e==2||e==4||e==6||e==7||e==9||e==11)n=31;if(e==3||e==5||e==8||e==10)n=30;if(e==1){n=28;t%4==0&&(n=29);t%100==0&&(n=28);t%400==0&&(n=29)}dIM=n}function frac(e){e-=Math.floor(e);e<0&&(e+=1);return e}function rnd(num,num2){with(Math)num=round(num*pow(10,num2))/pow(10,num2);return num}var dat,JD,UT,offset,dIM,RA,EOT,year,month,day,hours,minutes,seconds,utYear,utMonth,utDay,utHours,lat,longit,offset,locOffset,monthName=new Array("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),dayName=new Array("Son","Mon","Tue","Wed","Thu","Fri","Sat"),UTRISE,UTSET,RISE,SETT,ABOVE,ok,OK,twilight,LatLong,locName,hRise1,hSet1,hSet1,hSet2,ns,ew,string1="ex1",maxElevString,transString,timezoneString,lengthDayString,deg=String.fromCharCode(176),strA=new Array(3),strB=new Array(25),strArray=new Array(31),text=new Array(8);