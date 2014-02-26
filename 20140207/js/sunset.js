/*
 JavaScript Sun Table Calculator
© 2001 Juergen Giesen
http://www.jgiesen.de
*/

var dat, JD, UT, offset, dIM, RA, EOT;
var year, month, day, hours, minutes, seconds, utYear, utMonth, utDay, utHours;
var lat, longit, offset, locOffset;
var monthName = new Array('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
var dayName = new Array('Son','Mon','Tue','Wed','Thu','Fri','Sat');

var UTRISE, UTSET, RISE, SETT, ABOVE, ok, OK, twilight, LatLong, locName;
var hRise1, hSet1, hSet1, hSet2, ns, ew;
var  string1="ex1", maxElevString,transString, timezoneString, lengthDayString;
var deg = String.fromCharCode(176);

var strA = new Array(3);
var strB = new Array(25);
var strArray = new Array(31);
var text = new Array(8);

function setMyMenu(d, m, y, h, min) {
 //if ((d>0) && (d<32))
  //document.myform.Tag.options[d-1].selected=true; else console.log("Error on day: "+d);
 //if ((m>0) && (m<13))
  //document.myform.Monat.options[m-1].selected=true; else console.log("Error on month: "+m);
	//document.myform.Jahr.options[y-1992].selected=true;
 //if ((h>=0) && (h<24))
	 //document.myform.Stunde.options[h].selected=true; else console.log("Error on hour: "+h);
	//document.myform.Minute.options[min].selected=true;
}

function currentTime() {

 var curDat = new Date();
 var curOffset = curDat.getTimezoneOffset();
 if (curOffset>=1380) curOffset=curOffset-1440;
 var dLocOffset=-curOffset/60;

 utHours = curDat.getHours() - dLocOffset;
 LatLong =  document.myform.location.options[document.myform.location.selectedIndex].value;
 if (document.myform.location.selectedIndex != 0) getLocation(LatLong);
 var nn=2*locOffset;
 if (locOffset>0) nn=nn-1;
 document.myform.timezone.options[Math.abs(nn)].selected=true;

 hours = utHours + locOffset;
 if (hours<0) hours+=24;
 if (hours>23) hours-=24;
 minutes = curDat.getMinutes();
 seconds = curDat.getSeconds();
 seconds = 0;
 year = curDat.getYear();
 if (year<1900) year=year+1900;
 month = curDat.getMonth()+1;
 day = curDat.getDate();
 utMinutes = minutes;
 utMonth = month;
 utYear = year;
 document.myform.localDateText.value = writeDateTime(hours,day,month);
 setMyMenu(day, month, year, hours, minutes);
}

function writeDateTime(h, d, m) {

 var wDat = new Date(year,m-1,d,hours,minutes,seconds);
 var dayOfWeek = dayName[wDat.getDay()];
 var str = " " + year+", " + dayOfWeek + " " + document.myform.Monat.options[m-1].value+" ";
 if (d<10) str=str+"0"+d; else str=str+d;
 str=str+" at ";
 if (h<10) str=str+"0"+h+":"; else  str=str+h+":";
 if (minutes<10) str=str+"0"+minutes; else str=str+minutes;
 return str;
}

function getDateTime() {

 //var LatLong = document.myform.location.options[document.myform.location.selectedIndex].value;
 var LatLong = getLocation("52.38/4.89/+1*Amsterdam");
 var n=LatLong.length;
 var latStr="", s=0, s1=0, star=0;
 for (var i=0; i<n; i++) {if (LatLong.charAt(i)=='/') {s=i;break}};
 for (var i=s+1; i<n; i++) {if (LatLong.charAt(i)=='/') {s1=i;break}};
 for (var i=s1+1; i<n; i++) {if (LatLong.charAt(i)=='*') {star=i;break}};
 locName = LatLong.substring(star+1,n);

 offset=-60*document.myform.timezone.options[document.myform.timezone.selectedIndex].value;
 locOffset = -offset/60;

 var dat = new Date();
 hours = dat.getHours();
 minutes = dat.getMinutes();
 seconds = dat.getSeconds();
 year = dat.getYear();
 if (year<1900) year=year+1900;
 month = dat.getMonth()+1;
 day = dat.getDate();

 var millis = dat.getTime();
 millis=millis-locOffset*3600*1000;
 var d = new Date(millis);
 utHours = d.getHours();
 utDay = d.getDate();
 utMonth = d.getMonth()+1;
 utYear = d.getYear();
 if (utYear<1900) utYear=utYear+1900;
 seconds = 0;
 UT = utHours + minutes/60 + seconds/3600;

 document.myform.localDateText.value = writeDateTime(hours,day,month);
 document.myform.UTimeText.value = writeDateTime(utHours,utDay,utMonth);
 setMyMenu(day, month, year, hours, minutes);

 EOT = eot(utDay,utMonth,utYear,UT);
 theJulDay();
}

function getTZone() {
 
 var hRISE1, hSET1, hRISE2, hSET2;
 var index=document.myform.timezone.selectedIndex;
 var delta = document.myform.timezone.options[index].value;

 var dat = new Date(Date.UTC(year,month-1,utDay,utHours,minutes,seconds));
 var millis = Date.UTC(year,month-1,day,utHours,minutes,seconds);
 var theOffset = dat.getTimezoneOffset();
 if (theOffset>=1380) theOffset=theOffset-1440;
 var dLocOffset=-theOffset/60;

 millis = millis - (dLocOffset -delta)*3600*1000;

 var dat1 = new Date(millis);
 hours = dat1.getHours();
 minutes = dat1.getMinutes();
 seconds = dat1.getSeconds();
 year = dat1.getYear();
 if (year<1900) year=year+1900;
 month = dat1.getMonth()+1;

 day = dat1.getDate();

 if (delta==0) {day=utDay;month=utMonth;}

 if (delta>0) {

 if ((month>utMonth) && (Number(utHours)+Number(delta)>=24))
 {month=utMonth; day=utDay+1; daysInMonth(month,year); if (day>dIM) {day=1;month=month+1;};};

 if ((month==utMonth) || (Number(utHours)+Number(delta)<24))
  {
  if (day>utDay+1) {day=day-1;};
  if (month>utMonth) {day=1};
  if (Number(utHours)+Number(delta)<24) {day=utDay; month=utMonth; year=utYear;}
  if ((year>utYear)  && (day>1)) {day=day-1;}
  }
 }

if (delta<0) {

 daysInMonth(month,year);
 if ((day==dIM) && (utDay==2)) {day=1;month=utMonth}; 

 if ((day!=dIM) && (utDay!=2)) {
 if (day<utDay-1) day=day+1;
  if (month<utMonth) {daysInMonth(month,year);day=dIM}
  if (Number(utHours)+Number(delta)>=0) {day=utDay; month=utMonth;}
  if (year<utYear)  {if (day<31) day=day+1;}
 }
}

 setMyMenu(day, month, year, hours, minutes);
 document.myform.localDateText.value = writeDateTime(hours,day,month);
 calcRiseSet();
}

function ut() {
 utHours = hours - locOffset;
}

function setUT() {
 
 seconds =0;
 var theDate = new Date(year,month-1,day,hours,minutes,seconds);
 var millis = theDate.getTime();

 millis=millis-locOffset*3600*1000;

 var d = new Date(millis);
 utHours = d.getHours();
 utMinutes = d.getMinutes();
 utDay = d.getDate();
 utMonth = d.getMonth()+1;
 utYear = d.getYear();
 if (utYear<1900) utYear=utYear+1900;

 UT = utHours + utMinutes/60 + seconds/3600;
 document.myform.localDateText.value = writeDateTime(hours,day,month);
 document.myform.UTimeText.value = writeDateTime(utHours,utDay,utMonth);
 EOT = eot(utDay,utMonth,utYear,UT);
}


function getLocation(LatLong) {

 var n=LatLong.length;
 var latStr="", s=0, s1=0, star=0;
 for (var i=0; i<n; i++) {if (LatLong.charAt(i)=='/') {s=i;break}};
 for (var i=s+1; i<n; i++) {if (LatLong.charAt(i)=='/') {s1=i;break}};
 for (var i=s1+1; i<n; i++) {if (LatLong.charAt(i)=='*') {star=i;break}};

 var latStr = LatLong.substring(0,s);
 //document.myform.latitude.value=Math.abs(latStr);
 lat = Number(latStr);
 var globalLat= lat;
 //if (lat>=0) {ns=" N";  document.myform.NorthSouth.options[0].selected=true;}
 //else {ns=" S"; document.myform.NorthSouth.options[1].selected=true;}

 var longStr = LatLong.substring(s+1,s1);
 //document.myform.longitude.value=Math.abs(longStr);
 longit = Number(longStr);
 var globalLong= longit;
 //if (longit>=0) {ew=" E"; document.myform.EastWest.options[0].selected=true;}
 //else {ew=" W"; document.myform.EastWest.options[1].selected=true;}
 var tzString = LatLong.substring(s1+1,star);
 locOffset=Number(tzString);
 var nn=2*locOffset;
 if (locOffset>0) nn=nn-1;
 //document.myform.timezone.options[Math.abs(nn)].selected=true;

 locName = LatLong.substring(star+1,n);

 var millis = Date.UTC(year,utMonth-1,utDay,utHours,minutes,seconds);
 millis=millis+locOffset*3600*1000;

 var d = new Date();
 var doffset=d.getTimezoneOffset();
 if (doffset>=1380) doffset=doffset-1440;
 var dlocOffset=-doffset/60;
 millis=millis-dlocOffset*3600*1000;
 d = new Date(millis);
 hours = d.getHours();
 day = d.getDate();
 month = d.getMonth()+1;
 year = d.getYear();
 if (year<1900) year=year+1900;

 setMyMenu(day, month, year, hours, minutes);

 if (locOffset>=0) timezoneString="GMT + " + locOffset;
 else  timezoneString="GMT  " + locOffset;

 UT = utHours + minutes/60 + seconds/3600;
 //document.myform.localDateText.value = writeDateTime(hours,day,month);

 calculate();
 calcRiseSet();
}

//-- quadratic approximation and iteration for rise and set --

function riseset(DATE, MONTH, YEAR, HOUR) {
 var K = Math.PI/180.0;
 var sh = Math.sin(K*twilight);
 var dec = declination(DATE, MONTH, YEAR, HOUR);
 var gha = computeGHA (DATE, MONTH, YEAR, HOUR)	
 Y0 = Math.sin(K*computeHeight(dec, lat, longit, gha)) - sh;
	
 dec = declination(DATE, MONTH, YEAR, HOUR+1);
 gha = computeGHA (DATE, MONTH, YEAR, HOUR+1.0)	
 yPlus = Math.sin(K*computeHeight(dec, lat, longit, gha)) - sh;
	
 dec = declination(DATE, MONTH, YEAR, HOUR-1);
 gha = computeGHA (DATE, MONTH, YEAR, HOUR-1.0)
 yMinus = Math.sin(K*computeHeight(dec, lat, longit, gha)) - sh;
	
	ABOVE = (yMinus>0.0);

	QUAD(yMinus,yPlus);
	switch (NZ)
	{
		case 0: break;
		case 1:		
			if (yMinus<0.0) {UTRISE = HOUR+zero1; RISE=true;}
			else {UTSET = HOUR+zero1; SETT=true;}
			break;		
		case 2:		
			if (YE<0.0) {UTRISE = HOUR+zero2; UTSET=HOUR+zero1;}
			else {UTRISE = HOUR+zero1; UTSET=HOUR+zero2;}
			RISE=true;
			SETT=true;
			break;					
	}		
}

function QUAD(yMinus,yPlus){
	NZ = 0;
	var XE;
	var A = 0.5*(yMinus+yPlus) - Y0;
	var B = 0.5*(yPlus-yMinus);
	var C = Y0;
	XE = -B/(2.0*A);
	YE = (A*XE+B)*XE + C;
	var DIS = B*B - 4.0*A*C;
	if (DIS>=0) {
		DX = 0.5*Math.sqrt(DIS)/Math.abs(A);
		zero1 = XE - DX;
		zero2 = XE + DX;
		if (Math.abs(zero1)<=1.0) NZ = NZ + 1;
		if (Math.abs(zero2)<=1.0) NZ = NZ + 1;
		if (zero1<-1.0) zero1=zero2;
	}
}

// -------- Astro functions --------------------------------------

function JulDay (d, m, y, u){
 if (y<1900) y=y+1900
 if (m<=2) {m=m+12; y=y-1}
 A = Math.floor(y/100);
 JD =  Math.floor(365.25*(y+4716)) + Math.floor(30.6001*(m+1)) + d - 13 -1524.5 + u/24.0;
 return JD
}

function computeHeight(dec, latitude, longit, gha) {	
 var K = Math.PI/180.0;
 var lat_K = latitude*K;
 var dec_K = dec*K;
 var x=Number(gha)+Number(longit);
 sinHeight =Math.sin(dec_K)*Math.sin(lat_K) + Math.cos(dec_K)*Math.cos(lat_K)*Math.cos(K*x);	
 return Math.asin(sinHeight)/K;
}

function declination(day, month, year, UT) {

 var K = Math.PI/180.0;
 var jd =JulDay(day, month, year, UT);
 var T = (jd - 2451545.0)/36525.0;
 var L0 = 280.46645 + (36000.76983 + 0.0003032*T)*T;
 var M = 357.52910 + (35999.05030 - (0.0001559*T + 0.00000048*T)*T)*T; 
 M = K*M;
 var C = (1.914600 - 0.004817*T - 0.000014*T*T)*Math.sin(M) + (0.019993 - 0.000101*T)*Math.sin(2*M) + 0.000290*Math.sin(3*M);
 var theta = L0 + C;
 var omega = 125.04 - 1934.136*T;
 var lambda = theta - 0.00569 - 0.00478*Math.sin(K*omega);
 var eps0 =  23.0 + 26.0/60.0 + 21.448/3600.0 - (46.8150*T + 0.00059*T*T - 0.001813*T*T*T)/3600;
 var eps = eps0 + 0.00256*Math.cos(K*omega);
 var declin = Math.sin(K*eps)*Math.sin(K*lambda);
 declin = Math.asin(declin)/K;
 RA =Math.atan2(Math.cos(K*eps)*Math.sin(K*lambda), Math.cos(K*lambda))/K;
 if (RA<0) RA=RA+360;
 return declin;
}

function computeGHA (T, M, J, STD) {	
		var K = Math.PI/180.0, N, X, XX, P;		
		 N = 365 * J + T + 31 * M - 46;		 
		 if (M<3) N = N + Math.floor((J-1)/4);
		else N = N - Math.floor(0.4*M + 2.3) + Math.floor(J/4);			 
		P = STD/24.0;
		X = (P + N - 7.22449E5) * 0.98564734 + 279.306;
		X = X * K;
		XX = -104.55 * Math.sin(X) - 429.266 * Math.cos(X) + 595.63 * Math.sin(2.0 * X) - 2.283 * Math.cos(2.0 * X);
		XX = XX + 4.6 * Math.sin(3.0 * X) + 18.7333 * Math.cos(3.0 * X);
		XX = XX - 13.2 * Math.sin(4.0 * X) - Math.cos(5.0 * X) - Math.sin(5.0 * X) / 3.0 + 0.5 * Math.sin(6.0 * X) + 0.231;
		XX = XX / 240.0 + 360.0 * (P + 0.5);
		if (XX > 360.0) XX = XX - 360.0;  
 	return XX;				
}

function computeAzimut(dec, latitude, longitude, GHA, hoehe) {	
		var K = Math.PI/180.0;
		var cosAz, Az,nenner;
		var lat_K = latitude*K;
		var hoehe_K = hoehe*K;		
		nenner = Math.cos(hoehe_K)*Math.cos(lat_K);	
			cosAz = (Math.sin(dec*K) - Math.sin(lat_K)*Math.sin(hoehe_K))/nenner;					
		Az = Math.PI/2.0 - Math.asin(cosAz);
		Az = Az/K;
		if (Math.sin(K*(Number(GHA)+Number(longitude))) <= 0) Az = Az;	
		else Az = 360.0 - Az;			
		return Az;
	}

function GM_Sidereal_Time (JD) {	
 var t_eph, ut, MJD0, MJD;				
 MJD = JD - 2400000.5;		
 MJD0 = Math.floor(MJD);
 ut = (MJD - MJD0)*24.0;		
 t_eph  = (MJD0-51544.5)/36525.0;			
 return  6.697374558 + 1.0027379093*ut + (8640184.812866 + (0.093104 - 0.0000062*t_eph)*t_eph)*t_eph/3600.0;		
}

function LM_Sidereal_Time (JD, longit) {
 var GMST = GM_Sidereal_Time(JD);			
 var LMST =  24.0*frac((GMST + longit/15.0)/24.0);
 //return HoursMinutes(LMST);
return HoursMinutesSeconds(LMST);
}

function sunL(T){
	var L = 280.46645 + 36000.76983*T + 0.0003032*T*T	
	L = L % 360		
	if (L<0) L = L + 360
	return L			
}

function deltaPSI(T){
	var K = Math.PI/180.0
	var deltaPsi, omega, LS, LM				
	LS = sunL(T)	
	LM = 218.3165 + 481267.8813*T		
	LM = LM % 360	
	if (LM<0) LM = LM + 360		
	omega = 125.04452 - 1934.136261*T + 0.0020708*T*T + T*T*T/450000
	deltaPsi = -17.2*Math.sin(K*omega) - 1.32*Math.sin(K*2*LS) - 0.23*Math.sin(K*2*LM) + 0.21*Math.sin(K*2*omega)
	deltaPsi = deltaPsi/3600.0		
	return deltaPsi	
}
	
function EPS(T) {
	var K = Math.PI/180.0
	var LS = sunL(T)
	var LM = 218.3165 + 481267.8813*T	
	var eps0 =  23.0 + 26.0/60.0 + 21.448/3600.0 - (46.8150*T + 0.00059*T*T - 0.001813*T*T*T)/3600
	var omega = 125.04452 - 1934.136261*T + 0.0020708*T*T + T*T*T/450000		
	var deltaEps = (9.20*Math.cos(K*omega) + 0.57*Math.cos(K*2*LS) + 0.10*Math.cos(K*2*LM) - 0.09*Math.cos(K*2*omega))/3600
	return eps0 + deltaEps
}
	
function eot(date, month, year, UT) {		
	var K = Math.PI/180.0
	var  T = (JulDay(date, month, year, UT) - 2451545.0) / 36525.0		
	var eps = EPS(T);
 var dummy = declination(date, month, year, UT); // to get RA
	var LS = sunL(T)
	var deltaPsi = deltaPSI(T)						
	var E = LS - 0.0057183 - RA + deltaPsi*Math.cos(K*eps)		
	if (E>5) E = E - 360.0		
	E = E*4; // deg. to min		
	E = Math.round(1000*E)/1000								
	return E		
}


function transit(d,m,y) {

 var dummy =  declination(d, m, y, utHours);
 var trans = 12.0 - eot(d,m,y,12)/60.0 - longit/15.0; 
 var GHA = computeGHA (d, m, y, trans);

 var declin = declination(d, m, y, trans)
 var elev = computeHeight(declin, lat, longit, GHA);
 maxElevString =  Math.round(10*elev)/10;
 if (Math.round(10*elev)/10==Math.round(elev)) maxElevString=maxElevString+".0";

 trans = trans+locOffset;
 transString = HoursMinutes(trans);
}

// ----------------------------------------------------------------------

function dayTable(what) {

 var str1, str2, str3;
 var s1, s2, s3, s4, s5, s6, s7;

 var str, declin, GHA, elev, azimuth, std;
 calcRiseSet();
 if (lat>=0) ns=" N"; else ns=" S";
 if (longit>=0) ew=" E"; else ew=" W";
 if (locOffset>=0) timezoneString="GMT + " + locOffset;
 else  timezoneString="GMT  " + locOffset;

 str = locName + ",  " + timezoneString + "\n";
 document.myform.area.value+=str;
 s1=str;
 str = "Latitude: " + Math.abs(lat) +deg + ns + ",  Longitude: " + Math.abs(longit)+deg + ew + "\n";
 document.myform.area.value+=str;
 s2=str;
 document.myform.area.value+="Local time: " + document.myform.localDateText.value + "\n";
 s3="Local time: " + document.myform.localDateText.value + "\n";

 str = "Sunrise: " +  document.myform.riseTime.value;
 str = str + "     Sunset: " + document.myform.setTime.value;
 str = str + "    Length of Day: " +  lengthDayString + " h";
 document.myform.area.value+=str + "\n";
 s4=str;
 str = "Equation of Time: " + EOT + " min = " + HoursMinutes(EOT) + " min";
 document.myform.area.value+=str + "\n";
 s5=str;
 var declin = declination(utDay, utMonth, utYear, UT);
 declin = Math.round(1000*declin)/1000;
 str = "Declination: " + declin + " deg = " + DegreesMinutes(declin);
 if (declin>=0) str=str+" N"; else str=str+" S";
 document.myform.area.value+=str + "\n";
 s6=str;
  transit(day,month,year);
 document.myform.area.value+= "Transit at: " + transString + "   Altitude: " + maxElevString + deg + "\n";
 s7= "Transit at: " + transString + "&nbsp;&nbsp;&nbsp;Altitude: " + maxElevString + deg;
 document.myform.area.value+="Table of Altitude, Azimuth, and GHA Angles:" + "\n";

 for (var i=-locOffset; i<=-locOffset+24; i++) {
 std = i+locOffset;

 declin = declination(utDay, utMonth, utYear, i);
 GHA = computeGHA (utDay, utMonth, utYear, i);
 elev = computeHeight(declin, lat, longit, GHA);
 azimuth = computeAzimut(declin, lat, longit, GHA, elev);

 azimuth = Math.round(10*azimuth)/10;
 if (azimuth<100) azimuth="  "+azimuth;
 if (azimuth<10) azimuth="   "+azimuth;
 if (Math.round(azimuth)==Math.round(10*azimuth)/10) azimuth=azimuth+".0";

 elev = Math.round(10*elev)/10;
 if (Math.round(elev)==Math.round(10*elev)/10) elev=elev+".0";

 if ((elev<10) && (elev>0)) elev=" "+elev;
 if (elev>0) elev="  "+elev;
 if ((elev>-10) && (elev<0)) elev=" "+elev;
 
 if (std<10) str1=" "+std; else str1=std;
 str1=str1+":00";
 str2 =elev;
 str3=azimuth;
 str4=Math.round(10*GHA)/10.0;
	if (Math.round(GHA)==Math.round(10*GHA)/10.0) str4=str4+".0"; 
 //var LHA = Number(GHA+longit);
 //var str5 = Math.round(10*LHA)/10.0;

 //if (what==1) document.myform.area.value+=str1 + "   " + str2 + "   " + str3 + "   " + str4 + "   " + str5 + "\n";
 if (what==1) document.myform.area.value+=str1 + "   " + str2 + "   " + str3 + "   " + str4 + "\n";
 else {
 strA = new Array(str1,str2+deg,str3+deg);
 strArray[std] = strA;
 text[0] = s1;
 text[1] = s2;
 text[2] = s3;
 text[3] = s4;
 text[4] = s5;
 text[5] = s6;
 text[6] = s7;
 }
}
 if (what==2) writeDayPage(strArray)
}


function table(what) {

 var str11="", str22="";
 var T;

 if (lat>=0) ns=" N"; else ns=" S";
 if (longit>=0) ew=" E"; else ew=" W";
 if (locOffset>=0) timezoneString="GMT + " + locOffset;
 else  timezoneString="GMT  " + locOffset;

 tDat = new Date();
 daysInMonth(month,year);
 var str = "Sunrise, Sunset, Length of Day, Transit, and Max. Altitude"+"\n";
 str = str +  "Latitude: " + Math.abs(lat) +deg + ns + ",  Longitude: " + Math.abs(longit) +deg + ew+ ",   " + locName + ",  " + timezoneString + "\n";
text[0] = "Latitude: " + Math.abs(lat) +deg + ns + ",  Longitude: " + Math.abs(longit) +deg + ew+ ",&nbsp;&nbsp;&nbsp;" + locName + ",  " + timezoneString;
 document.myform.area.value+=str;

var url = this.location.toString();

 if ((url.lastIndexOf(String.fromCharCode(106,103,105,101,115))==-1)  && (url.lastIndexOf(string1)==-1)) {document.myform.area.value+=ok+"\n";dIM=1};

 for (var i=1; i<=dIM; i++) {

   JD=JulDay (i, month, year, 12)
   EOT = eot(i, month, year, 12)
   day=i
   RISE = false;
   SETT = false;
   twilight = -0.8333;
			
   for (var j=-locOffset; j<-locOffset+24; j++) { // i=GMT							
      riseset(i, month, year, j)								
      if (RISE && SETT)  break;			
   }

 if (RISE || SETT) {
 if (RISE) {UTRISE=UTRISE + locOffset;if (UTRISE>=24) UTRISE -=24;if (UTRISE<0) UTRISE +=24;}
 if (SETT) {UTSET=UTSET + locOffset;if (UTSET>=24) UTSET -=24;if (UTSET<0) UTSET +=24;}
 lengthOfDay=UTSET-UTRISE;
 }

 T = UTRISE-locOffset+hours-utHours;
 if (T<0) T +=24;
 if (T>=24) T -=24;
 if (RISE) str1=HoursMinutes(T);
 else {if (ABOVE) str1="visible"; else str1="--.--";}		

 T = UTSET-locOffset+hours-utHours;
 if (T<0) T +=24;
 if (T>=24) T -=24;
 if (SETT) str2=HoursMinutes(T);
 else {if (ABOVE) {str2="visible"; lengthOfDay=24} else {str2="--.--";lengthOfDay=0}}

   RISE = false;
   SETT = false;
   twilight = -6.0;
			
   for (var j=-locOffset; j<-locOffset+24; j++) { // i=GMT							
      riseset(i, month, year, j)								
      if (RISE && SETT)  break;			
   }

 if (RISE || SETT) {
 if (RISE) {UTRISE=UTRISE + locOffset;if (UTRISE>=24) UTRISE -=24;if (UTRISE<0) UTRISE +=24;}
 if (SETT) {UTSET=UTSET + locOffset;if (UTSET>=24) UTSET -=24;if (UTSET<0) UTSET +=24;}
 }

 T = UTRISE-locOffset+hours-utHours;
 if (T<0) T +=24;
 if (T>=24) T -=24;
 if (RISE) str11=HoursMinutes(T);
 else {if (ABOVE) str11="visible"; else str11="--.--";}		
	
 T = UTSET-locOffset+hours-utHours;
 if (T<0) T +=24;
 if (T>=24) T -=24;
 if (SETT) str22=HoursMinutes(T);
 else {if (ABOVE) str22="visible"; else str22="--.--";}		

 st1= ""+Math.round(100*lengthOfDay)/100;
 st2= ""+Math.round(10*lengthOfDay)/10;
 a= Math.round(100*lengthOfDay)/100;

 if (Math.round(lengthOfDay)==a) str=a+".00 h"; 
 else {if (st1.length==st2.length) str= a+"0 h"; else str=a+ " h"};
 
 if (i<10) s="0"; else s="";

 transit(day,month,year);
 str3 = transString;
 str4 = maxElevString+deg;

 var dummy =  declination(day,month,year, utHours);// fuer RA
 EOT = eot(day, month, year, 12+locOffset);
 str5 = DegreesMinutes(EOT);

 if (what==1)
 document.myform.area.value+=year + " " +monthName[month-1]+" "+s+day+"    "+str1 + "    " + str2 + "   " + str +"   "+str3+"  "+str4+"\n"

 strB = new Array(year,monthName[month-1],day,str11,str1,str2,str22,str,str3,str4,str5);
 strArray[i] = strB;

} // for day

 if (what==2) writeMonthPage(strArray,dIM);
}


function calcRiseSet () {
 offset=-60;
 locOffset = -offset/60;

 twilight = -0.8333;
 RISE = false;
 SETT = false;
 ok = OK.charAt(1)+OK.charAt(3)+OK.charAt(2)+OK.charAt(0);
 
 //document.myform.riseTime.value=ok;
 //document.myform.setTime.value=ok;
 //document.myform.tRiseTime.value=ok;
 //document.myform.tSetTime.value=ok;
		
 for (var i=-locOffset; i<-locOffset+24; i++) { // i=GMT							
  riseset(day, month,year, i)										
  if (RISE && SETT)  break;			
 }	

 if (RISE || SETT) {
 if (RISE) {
 hRise1=UTRISE + locOffset;
 if (hRise1>24) hRise1-=24;
 if (hRise1<0)   hRise1+=24;
 var declin = declination(utDay, utMonth, utYear, UTRISE);
 var GHA = computeGHA (utDay, utMonth, utYear, UTRISE);
 var elev = computeHeight(declin, lat, longit, GHA);
 var azimuth = computeAzimut(declin, lat, longit, GHA, elev);
 azimuth = Math.round(10*azimuth)/10;
 document.myform.azimRise.value=" "+azimuth;
 }

 if (SETT) {
 hSet1=UTSET + locOffset;
 if (hSet1>24) hSet1-=24;
 if (hSet1<0) hSet1+=24;
 var declin = declination(utDay, utMonth, utYear, UTSET);
 var GHA = computeGHA (utDay, utMonth, utYear, UTSET);
 var elev = computeHeight(declin, lat, longit, GHA);
 var azimuth = computeAzimut(declin, lat, longit, GHA, elev);
 azimuth = Math.round(10*azimuth)/10;
 //document.myform.azimSet.value=" "+azimuth;
 }
 lengthOfDay=hSet1-hRise1;
} // if (RISE || SETT)

var url = this.location.toString();
if ((url.lastIndexOf(String.fromCharCode(106,103,105,101,115))!=-1)  || (url.lastIndexOf(string1)!=-1))   {
 if (RISE) str=HoursMinutes(hRise1); else {if (ABOVE) str="visible"; else str="--.--";}			
 document.myform.riseTime.value=" "+str;	

 if (SETT) str=HoursMinutes(hSet1);
 else {if (ABOVE) {str="visible"; lengthOfDay=24;}  else {str="--.--";lengthOfDay=0}}			
 //document.myform.setTime.value=" "+str;
}
 lengthOfDay=Math.round(100*lengthOfDay)/100;
 //document.myform.text4.value=" "+lengthOfDay;
 lengthDayString=lengthOfDay;

 twilight = document.myform.twilight.options[document.myform.twilight.selectedIndex].value;
 RISE = false;
 SETT = false;
		
 for (var i=-locOffset; i<-locOffset+24; i++) { // i=GMT							
  riseset(day, month,year, i)										
  if (RISE && SETT)  break;			
 }

 if (RISE || SETT) {
 if (RISE) {hRise2=UTRISE + locOffset;if (hRise2>24) hRise2-=24;if (hRise2<0) hRise2+=24;}
 if (SETT) {hSet2=UTSET + locOffset;if (hSet2>24) hSet2-=24;if (hSet2<0) hSet2+=24;}
 lengthOfDay=hSet2-hRise2;
 }

 if ((url.lastIndexOf(String.fromCharCode(106,103,105,101,115))!=-1)  || (url.lastIndexOf(string1)!=-1))   {
 if (RISE) str=HoursMinutes(hRise2); else {if (ABOVE) str="visible"; else str="--.--";}		
 document.myform.tRiseTime.value=" "+str;

 if (SETT) str=HoursMinutes(hSet2);
 else {if (ABOVE) {str="visible"; lengthOfDay=24;}  else {str="--.--";lengthOfDay=0}}			

 document.myform.tSetTime.value=" "+str;
 }

 transit(day, month, year);
 document.myform.transText.value=" "+transString;
 document.myform.maxAltText.value=" "+maxElevString;
 
}

function calculate() {
 var globalLat, globalLong;
 var str = " " + Math.round(1000*declination(utDay,utMonth,utYear,UT))/1000;
 //document.myform.declinText.value=str;
 OK="ODMEDO";

// lat = document.myform.latitude.value;

// longit = document.myform.longitude.value;

 lat= globalLat;
 longit= globalLong;

 //if (document.myform.EastWest.selectedIndex==1) longit=-longit;
 //if (document.myform.NorthSouth.selectedIndex==1) lat=-lat;

 var declin = declination(utDay, utMonth, utYear, UT);
 var GHA = computeGHA (utDay, utMonth, utYear, UT);
 var elev = computeHeight(declin, lat, longit, GHA);
 elev =  Math.round(10*elev)/10;
 //document.myform.elevation.value=" "+elev;

 EOT = eot(utDay,utMonth,utYear,UT);
 str = " " + EOT;
 //document.myform.eotText.value=str;

 var azimuth = computeAzimut(declin, lat, longit, GHA, elev);
 //document.myform.azim.value = " "+Math.round(10*azimuth)/10;
 //document.myform.LMST.value=LM_Sidereal_Time(JulDay (day, month, year, UT),longit);

 GHA =  Math.round(10*GHA)/10;
 //document.myform.ghaText.value=" "+GHA;
 var ra = RA/15.0;
 if (ra>=24) ra=ra-24;
 //document.myform.raText.value=HoursMinutesSeconds(ra);
}


function theDay() {
 daysInMonth(month,year);
 var index=document.myform.Tag.selectedIndex;
 if (index<dIM) day= index + 1;
 else {day= dIM; document.myform.Tag.options[dIM-1].selected=true};
 setUT();
 //calculate();
 //calcRiseSet ();
}

function theMonth() {
 month = Number(document.myform.Monat.selectedIndex)+Number(1);
 theDay();
 //calculate();
 //calcRiseSet ();
}

function theYear() {
 year = document.myform.Jahr.selectedIndex + 1992;
 theDay();
 //calculate();
 //calcRiseSet ();
}

function theHour() {
 hours = document.myform.Stunde.selectedIndex;
}


function theMinute() {
 minutes = document.myform.Minute.selectedIndex;
 seconds = 0;
}

function theJulDay() {
 jd = JulDay (utDay, utMonth, utYear, UT)
	document.myform.JulDayText.value=" " + Math.round(10000*jd)/10000	
}

function theDateTime() {

 theDay();
 theMonth();
 theYear();
 theHour();
 theMinute();
 setUT();
 theJulDay();
 calculate();
 calcRiseSet();
}

function getLatitude() {
 locName="User Input";
 str=document.myform.latitude.value;
 var n=str.length;
 for (var i=0; i<n; i++) {
 c=str.charAt(i);
 if ((c!='0') && (c!='1')  && (c!='2') && (c!='3') && (c!='4') && (c!='5') && (c!='6') && (c!='7') && (c!='8') && (c!='9')  && (c!='+') && (c!='.')) {console.log("Error on latitude value !"+"\n"+"Enter positive decimal degree value, e.g.:"+"\n"+" 52.34, and select North or South.");document.myform.latitude.value=0; break;};}
 if (Math.abs(Number(str))>90) {console.log("Latitude must less or equal to 90 degrees !"); document.myform.latitude.value=0; };
 lat=  Number(str);
 if (document.myform.NorthSouth.selectedIndex==0) lat=Math.abs(lat);
 else lat=-Math.abs(lat);
 if (lat>=0) ns=" N"; else ns=" S";
 document.myform.location.options[0].selected=true;
}

function getLongitude() {
 str=document.myform.longitude.value;
 var n=str.length;
 for (var i=0; i<n; i++) {
 c=str.charAt(i);
 if ((c!='0') && (c!='1')  && (c!='2') && (c!='3') && (c!='4') && (c!='5') && (c!='6') && (c!='7') && (c!='8') && (c!='9')  && (c!='+') && (c!='.')) {console.log("Error on longitude value !"+"\n"+"Enter positive decimal degree value, e.g.:"+"\n"+" 8.34 and select East or West.");document.myform.longitude.value=0; break;};}
 if (Math.abs(Number(str))>180) {console.log("Longitude must less or equal to 180 degrees !"); document.myform.longitude.value=0; };
 longit=  Number(str);
 if (document.myform.EastWest.selectedIndex==0) longit=Math.abs(longit);
 else  longit=-Math.abs(longit);
 if (longit>=0) ew=" E"; else ew=" W";
 document.myform.location.options[0].selected=true;
 calculate();
 calcRiseSet();
}
function getTwilight() {
twilight = document.myform.twilight.options[document.myform.twilight.selectedIndex].value;
theDateTime();
}

function clearTable() {
 document.myform.area.value=" ";
}


function NS() {
 if (document.myform.NorthSouth.selectedIndex==0) lat=Math.abs(lat);
 else  lat=-Math.abs(lat);
 locName="User Input";
 document.myform.location.options[0].selected=true;
 calcRiseSet();
}

function EW() {
 if (document.myform.EastWest.selectedIndex==0) longit=Math.abs(longit);
 else  longit=-Math.abs(longit);
 locName="User Input";
 document.myform.location.options[0].selected=true;
 locName="User Input";
 calculate();
 calcRiseSet();
}

function yearTable() {

 var strM1 = new Array(13); var strM2 = new Array(13);
 var T;

 if (lat>=0) ns=" N"; else ns=" S";
 if (longit>=0) ew=" E"; else ew=" W";
 if (locOffset>=0) timezoneString="GMT + " + locOffset;
 else  timezoneString="GMT  " + locOffset;

 tDat = new Date();

 text[0] = year + ",&nbsp;&nbsp;&nbsp;" + "Latitude: " + Math.abs(lat) +deg + ns + ",  Longitude: " + Math.abs(longit) +deg + ew+ ",&nbsp;&nbsp;&nbsp;" + locName + ",  " + timezoneString;
 twilight = -0.8333;

  for (var i=1; i<=31; i++) {

 for (var mm=1; mm<=12; mm++) {

 daysInMonth(mm,year);

 RISE = false;
 SETT = false;
			
 for (var j=-locOffset; j<-locOffset+24; j++) { // i=GMT							
   riseset(i, mm, year, j)								
   if (RISE && SETT)  break;			
 }

 if (RISE || SETT) {
 if (RISE) {UTRISE=UTRISE + locOffset;if (UTRISE>=24) UTRISE -=24;if (UTRISE<0) UTRISE +=24;}
 if (SETT) {UTSET=UTSET + locOffset;if (UTSET>=24) UTSET -=24;if (UTSET<0) UTSET +=24;}
 }


 T = UTRISE-locOffset+hours-utHours;
 if (T<0) T +=24;
 if (T>=24) T -=24;

 if (RISE) str1=HoursMinutes(T);
 else {if (ABOVE) str1="vis."; else str1="--.--";}		

 T = UTSET-locOffset+hours-utHours;
 if (T<0) T +=24;
 if (T>=24) T -=24;

 if (SETT) str2=HoursMinutes(T);
 else {if (ABOVE) str2="vis."; else str2="--.--";}

 if (i<=dIM) {strM1[mm]=str1; strM2[mm]=str2;} else {strM1[mm]=""; strM2[mm]="";}

 var url = this.location.toString();
 if ((url.lastIndexOf(String.fromCharCode(106,103,105,101,115))==-1)  && (url.lastIndexOf(string1)==-1)) strM2[mm]=ok;

 } // for month


 strB = new  Array(i,strM1[1],strM2[1],strM1[2],strM2[2],strM1[3],strM2[3],strM1[4],strM2[4],strM1[5],strM2[5],strM1[6],strM2[6],strM1[7],strM2[7],strM1[8],strM2[8],strM1[9],strM2[9],strM1[10],strM2[10],strM1[11],strM2[11],strM1[12],strM2[12]);

strArray[i] = strB;

} // for day

 writeYearPage(strArray);
}

// -------------- Utility functions --------------------------------------

function HoursMinutes(time) {
 var t=time;
 time=Math.abs(time);
 var min=Math.round(60*(time-Math.floor(time)));
 var str;
 str = Math.floor(time);
 if (str<10) str="0"+str;
 if (min>=10) str=str+":"+min;
 else  str=str+":0"+min;
 if (min==60) str=Math.floor(time+1)+":00";
 if (t<0) return "-"+str; else return str;
}

function HoursMinutesSeconds(time) {
 var h = Math.floor(time);
 var min = Math.floor(60.0*frac(time));
 var secs = Math.round(60.0*(60.0*frac(time)-min));
 var str;
 if (min>=10) str=h+":"+min;
 else  str=h+":0"+min;
 //if (min==60) str=(h+1)+":00";
 if (secs<10) str = str + ":0"+secs;
 else str = str + ":"+secs;
 return " " + str;
}

function DegreesMinutes(time) {
 var t=time;
 time=Math.abs(time);
 var min=Math.round(60*(time-Math.floor(time)));
 var str;
 //if (min>=10) str=Math.floor(time)+" deg "+min + " min";
 //else  str=Math.floor(time)+" deg 0"+min + " min";
 if (min>=10) str=Math.floor(time)+deg+" "+min + "'";
 else  str=Math.floor(time)+deg+" 0"+min + "'";
 if (min==60) str=Math.floor(time+1)+":00" + "'";
 if (t<0) return "-"+str; else return str;
}

function daysInMonth(m, y) {
	var n=31
	m=m-1
	if ((m==0) || (m==2) || (m==4) || (m==6) || (m==7) || (m==9) || (m==11))  n=31
	if ((m==3) || (m==5) || (m==8) || (m==10))  n=30;
	if (m==1) {
		n=28;
		if ((y % 4) == 0) n=29
		if ((y % 100) == 0) n=28
		if ((y % 400) == 0) n=29
		}
	dIM=n;			
}


function frac(X) {
 X = X - Math.floor(X);
 if (X<0) X = X + 1.0;
 return X;		
}

function rnd(num,num2) {
  with (Math) {
    num = round(num*pow(10,num2)) / pow(10,num2)
  }
  return num
}


