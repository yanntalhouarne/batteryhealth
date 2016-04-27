#include "LiquidCrystal/LiquidCrystal.h"
#include "elapsedMillis/elapsedMillis.h"
#include "ThingSpeak/ThingSpeak.h"
#include <stdio.h>
#ifdef SPARK
	#include "ThingSpeak/ThingSpeak.h"
#else
	#include "ThingSpeak.h"
#endif

// On Particle: 0 - 4095 maps to 0 - 3.3 volts
#define VOLTAGE_MAX 3.3
#define VOLTAGE_MAXCOUNTS 4095.0

////////////      Counter  Globals          /////////////////
elapsedMillis timeElapsed1;
elapsedMillis timeElapsed2;
elapsedMillis timeElapsed3;
elapsedMillis timeElapsed4;
elapsedMillis timeElapsed5;
////////////      LCD  Globals          /////////////////

LiquidCrystal lcd(4, 2, 1, 3, 5, 0); // assign pins to LCD

byte z = 0b00000;           // special character definition (bars)
byte o = 0b11111;           //
byte bar0[7] = {z,z,z,z,o}; //
byte bar1[7] = {z,z,z,o,o}; //
byte bar2[7] = {z,z,o,o,o}; //
byte bar3[7] = {z,o,o,o,o}; //
byte bar4[7] = {o,o,o,o,o}; //

////////////      WiFi  Globals          /////////////////
int rssi = 0;
 TCPClient client;

//// Measurement globals ////

  double vv[1024], tt[1024];
  int nn = 0, ii;
  float soc_old, v, voc, capacity = 0, q, flag, current, soc_new;
  int mm = 0;
  int RawCurrent[500],RawVoltage[500];
  int Rest1 = 60000,Rest2 = 3600000, StoreTime = 600000;
  float voc_new,voc_old;
  int capPointer1 = 0;
  double capacityArray[1000];
  int capPointer=0; 
  unsigned long myChannelNumber = 86596;
  const char * myWriteAPIKey = "BXMUPPBX7V13ZRA5";

////////////////////////////////////////////////////////
////  function name:   getsoc()                     ////
////  paramater type: float                         ////
////  return type: float                            ////
////  summary: Lookup table to obtain SoC from VoC) ////
////////////////////////////////////////////////////////
float getsoc(float voc){
	int n=0;
	float socdata[100]={0,0.0101001852895826,
    0.0202003705791651,0.0303005558687477,
    0.0404007411583303,0.0505009264479129,
    0.0606011117374954,0.0707012970270780,
    0.0808014823166606,0.0909016676062431,
    0.101001852895826,0.111102038185408,
    0.121202223474991,0.131302408764573,
    0.141402594054156,0.151502779343739,
    0.161602964633321,0.171703149922904,
    0.181803335212486,0.191903520502069,
    0.202003705791651,0.212103891081234,
    0.222204076370817,0.232304261660399,
    0.242404446949982,0.252504632239564,
    0.262604817529147,0.272705002818729,
    0.282805188108312,0.292905373397895,
    0.303005558687477,0.313105743977060,
    0.323205929266642,0.333306114556225,
    0.343406299845807,0.353506485135390,
    0.363606670424973,0.373706855714555,
    0.383807041004138,0.393907226293720,
    0.404007411583303,0.414107596872885,
    0.424207782162468,0.434307967452051,
    0.444408152741633,0.454508338031216,
    0.464608523320798,0.474708708610381,
    0.484808893899963,0.494909079189546,
    0.505009264479129,0.515109449768711,
    0.525209635058294,0.535309820347876,
    0.545410005637459,0.555510190927041,
    0.565610376216624,0.575710561506207,
    0.585810746795789,0.595910932085372,
    0.606011117374954,0.616111302664537,
    0.626211487954119,0.636311673243702,
    0.646411858533285,0.656512043822867,
    0.666612229112450,0.676712414402032,
    0.686812599691615,0.696912784981197,
    0.707012970270780,0.717113155560363,
    0.727213340849945,0.737313526139528,
    0.747413711429110,0.757513896718693,
    0.767614082008275,0.777714267297858,
    0.787814452587441,0.797914637877023,
    0.808014823166606,0.818115008456188,
    0.828215193745771,0.838315379035353,
    0.848415564324936,0.858515749614519,
    0.868615934904101,0.878716120193684,
    0.888816305483266,0.898916490772849,
    0.909016676062431,0.919116861352014,
    0.929217046641597,0.939317231931179,
    0.949417417220762,0.959517602510344,
    0.969617787799927,0.979717973089509,
    0.989818158379092,0.999918343668675};
    
	float vocdata[100]={22.2644785478548,
    22.4817174787864,22.6241370260014,
    22.7256492977894,22.8155350068421,
    22.9064951804524,23.0000023194507,
    23.0294916092603,23.0778816234939,
    23.1390387089890,23.1965626498040,
    23.2489134032897,23.2983559726528,
    23.3460298311019,23.3929090361848,
    23.4458952230590,23.5040368783206,
    23.5301188493632,23.5472168140619,
    23.5882263523390,23.6286666351280,
    23.6683659526309,23.7128308691334,
    23.7599441994887,23.7817097034770,
    23.8149383549895,23.8524085092930,
    23.8890738577804,23.9249026124686,
    23.9601081052575,23.9953867502336,
    24.0316731863818,24.0765889919828,
    24.1109938946420,24.1293709725489,
    24.1635791845702,24.1983864377043,
    24.2336238740367,24.2675650819195,
    24.3012847564877,24.3345248333331,
    24.3679433127811,24.4030364629443,
    24.4369959106159,24.4700833981975,
    24.5029960578957,24.5355126690921,
    24.5683123774126,24.6016598817068,
    24.6353411224374,24.6679298425516,
    24.7001061644399,24.7320869545049,
    24.7639528103238,24.7966865509395,
    24.8296669599982,24.8619795620849,
    24.8937110473039,24.9250814035122,
    24.9567386602438,24.9888634611142,
    25.0203399515818,25.0515823992883,
    25.0824512878667,25.1132177035289,
    25.1440758507019,25.1752072046624,
    25.2066514672237,25.2372673764743,
    25.2676237365192,25.2979478548321,
    25.3281949118856,25.3595624367075,
    25.3908075017186,25.4208417531836,
    25.4510349786246,25.4810428165172,
    25.5109723899603,25.5410009102072,
    25.5718528616999,25.6020936158728,
    25.6313317234860,25.6612723776144,
    25.6906940093807,25.7203946409677,
    25.7509983038171,25.7812022719959,
    25.8147209651078,25.8481724401536,
    25.8904009187668,25.9343318964089,
    25.9787549684130,26.0239415322818,
    26.0751377146438,26.1305531679651,
    26.1919572043993,26.2614318138621,
    26.3369100903683,26.4335607670413,
    26.5401697290705};
    while (n<99)
    {
	    if (voc<=vocdata[n])
        {
		    break;
	    }
	n++;
    }	
return socdata[n];	
}
/////////////     End of getsoc()       /////////////////



///////////////////////////////////////////////////////////////
////  function name:   Wifi_Strength()                     ////
////  paramater type: void                                 ////
////  return type: void                                    ////
////  summary: Dipplays bars depending on WiFI connetivity ////
///////////////////////////////////////////////////////////////
void WiFi_Strength()
{
  // WiFi strength display

      //strength
            lcd.setCursor(12, 1);
            if ( rssi < -80) 
            {
                lcd.write(0);
            }
            if ( rssi < -60) 
            {
                lcd.setCursor(12, 0);
                lcd.write(0);
                lcd.setCursor(13, 0);
                lcd.write(1);
            }
            else if ( rssi < -50 )
            { 
                lcd.setCursor(12, 0);
                lcd.write(0);
                lcd.setCursor(13, 0);
                lcd.write(1);
                lcd.setCursor(14, 0);
                lcd.write(2);
                
            }
           else if ( rssi < -20 )
            { 
                lcd.setCursor(12, 0);
                lcd.write(0);
                lcd.setCursor(13, 0);
                lcd.write(1);
                lcd.setCursor(14, 0);
                lcd.write(2);
                lcd.setCursor(15, 0);
                lcd.write(3);
            }
           else if ( rssi > -19  &&  rssi < 0 )
            { 
                lcd.setCursor(12, 0);
                lcd.write(0);
                lcd.setCursor(13, 0);
                lcd.write(1);
                lcd.setCursor(14, 0);
                lcd.write(2);
                lcd.setCursor(15, 0);
                lcd.write(3);
                lcd.setCursor(16, 0);
                lcd.write(4);
            }
            else if (rssi >0) 
            {
              lcd.setCursor(12, 0);
              lcd.print("OFF  ");
            }
}
/////////////    en of Wifi_Strength()       /////////////////


///////////////////////////////////////////////////////////////////////
////  function name:   bars_characters()                           ////
////  paramater type: void                                         ////
////  return type: void                                            ////
////  summary: Creates the bar characters for the WiFi onnectivity ////
///////////////////////////////////////////////////////////////////////
void bars_characters()
{
    lcd.createChar(0, bar0);
    lcd.createChar(1, bar1);
    lcd.createChar(2, bar2);
    lcd.createChar(3, bar3);
    lcd.createChar(4, bar4);
}
/////////////    end of bars_charcters       //////////////////////////


///////////////////////////////////////////////////////////////////////
////  function name:   corr()                                      ////
////  paramater type: double *x, double *y,int n                   ////
////  return type: double                                          ////
////  summary:                                                     ////
///////////////////////////////////////////////////////////////////////
double corr(double *x, double *y,int n)
{
	double	EX = 0;
	double EY = 0;
	double EXY = 0;
	double EXsq = 0;
	double EYsq = 0;
	for (int i = 0; i < n; i++)
	{
		EX = EX + x[i];
		EY = EY + y[i];
		EXY = EXY + x[i]*y[i];
		EXsq = EXsq + x[i] * x[i];
		EYsq = EYsq + y[i] * y[i];
	}
	return (n*EXY - EX*EY) / sqrt((n*EXsq - EX*EX)*(n*EYsq - EY*EY));
}
////////////////////////////    end of corr()       ///////////////////


///////////////////////////////////////////////////////////////////////
////  function name:   OCVs()                                      ////
////  paramater type: double* t, double* v, int n                  ////
////  return type: double                                          ////
////  summary:                                                     ////
///////////////////////////////////////////////////////////////////////
double OCVss(double* t, double* v, int n)
{
	double vmax = 37;
	double vmin = 18;
	double* T = new double[n];
	double* V = new double[n];
	for (int i = 0; i < n; i++)
	{
		T[i] = t[i];
		V[i] = v[i];
	}
	int initGuessNum = 101;
	double *guess = new double[initGuessNum];
	double *corrArray = new double[initGuessNum];
	double inc = (vmax - vmin) /(double)(initGuessNum - 1);
	double* vp = new double[n];
	
	for (int i = 0; i < initGuessNum; i++)
	{
		guess[i] = vmin + inc*i;
	//	cout << guess[i] << endl;
		int j;
		for (j = 0; j < n; j++)
		{
			vp[j] = -V[j] + guess[i];
			if (vp[j] <= 0)
			{
				break;
			}
			else
			{
				vp[j] = log(vp[j]);
			}
		}
		if (j != n)
		{
			if (j < 10)
			{
				corrArray[i] = 0;
			}
			else
			{
				double* Ttc = new double[j];
				double* vptc = new double[j];
				for (int k = 0; k < j; k++)
				{
					Ttc[k] = T[k];
					vptc[k] = vp[k];
				}
				corrArray[i] = corr(vptc, Ttc, j);
                delete Ttc;
                delete vptc;
			}
		}
		else
		{
			corrArray[i] = corr(vp, T, n);
		}
	}
	int minIndex = 0;
	double minCorr = corrArray[0];
	//cout << "corrArray" << endl;
	for (int i = 1; i < initGuessNum; i++)
	{
		//cout << corrArray[i] << endl;
		if (corrArray[i] < minCorr)
		{
			minCorr = corrArray[i];
			minIndex = i;
		}
	}
	vmax = guess[minIndex + 1];
	vmin = guess[minIndex - 1];
	double retVal=guess[minIndex];
	delete T;
	delete V;
	delete guess;
	delete corrArray;
	delete vp;
	return retVal;
}
////////////////////////   end of OCVss()       ///////////////////////

void setup() 
{
  ThingSpeak.begin(client); // call ThingSpeak
 
  Particle.variable("RSSI",&rssi,INT); // return rssi variable for connectivity information
  
  lcd.begin(16, 2); // turn on LCD
  
  lcd.print("Initializing..."); //Initialize LCD
  delay(1000);
  lcd.clear(); 
  lcd.setCursor(0, 0); //set curor to left

  bars_characters(); //create the bar characters
}

   
void loop() {
    //get the rssi varible which indicates the WiFi connectivity status    
     rssi = WiFi.RSSI();
     Particle.publish("rssi",String(rssi),60,PRIVATE);
    
     WiFi_Strength(); //update the WiFi connectivity as bars on CD
    
     /*RawVoltage[mm] = analogRead(A0)*3.3*11/4096;
     RawCurrent[mm] = analogRead(A1)*0.082352 - 164.786;
     v = 0;
     current = 0;
     mm = (mm + 1) % 500;
     for(int z = 0; z < 500 ; z++)
     {
         v = v + RawVoltage[z];
        current = current + RawCurrent[z];
     }
     v = v/500.0;
     current=current/500.0;
     */
     v = analogRead(A0)*3.3*11/4096;
     current = analogRead(A1)*0.082352 - 164.786;
     current = -current;
    
    if(current < -5)       //charging state
    {
        //lcd.setCursor(0, 0);
        //lcd.print("Charging");
        soc_new = soc_old = voc_old = voc_new = flag = q = 0;
        timeElapsed1 = timeElapsed2 = timeElapsed3 = timeElapsed4 = 0;
        for (ii = 0; ii <= nn; ii++)
            {
                vv[ii] = 0;
                tt[ii] = 0;
            }
            nn = 0;
    }
    else
    {
    if (current > 20)       //working state
    {
     //lcd.setCursor(0, 0);
     //lcd.print("Working");
     
        if (timeElapsed1 > Rest1 && timeElapsed1 < Rest2)
        {
            if (timeElapsed1 > Rest1 + StoreTime)
            {
                if(voc_new >= 0.01)
                {
                    voc_old = voc_new;
                }
                voc_new = OCVss(tt,vv,nn);
                if(voc_old > 0.01)
                {
                    soc_old = getsoc(voc_old);
                    soc_new = getsoc(voc_new);
                    capacity = q/(soc_old - soc_new);
                    capacityArray[capPointer++] = capacity;
                    ThingSpeak.writeField(myChannelNumber, 3, capacity, myWriteAPIKey);
                }
                q = 0;
                
            }
            for (ii = 0; ii <= nn; ii++)
            {
                vv[ii] = 0;
                tt[ii] = 0;
            }
            nn = 0;
        }
        timeElapsed1 = 0;
        flag = 0;
    }
    else
    {
     //lcd.setCursor(0, 0);
     //lcd.print("Resting");
     
    }
    
    if (timeElapsed1 > Rest1 && !flag)
    {
        if(timeElapsed5 > 5000)
        {
            timeElapsed5 = 0;
            vv[nn] = v;
            tt[nn] = timeElapsed1/1000;
            nn++;
        }
        if (timeElapsed1 > Rest2)
        {
            if(voc_new >= 0.01)
            {
                voc_old = voc_new;
            }
            voc_new = v;
            if(voc_old > 0.01)
            {
                soc_old = getsoc(voc_old);
                soc_new = getsoc(voc_new);
                capacity = q/(soc_old - soc_new);
                capacityArray[capPointer++] = capacity;
                ThingSpeak.writeField(myChannelNumber, 3, capacity, myWriteAPIKey);
            }
            q = 0;
            flag = 1;
        }
    }
    q = q + current * timeElapsed2*1.0/1000;
    timeElapsed2 = 0;
    

         
    /*     lcd.setCursor(0,0);
         lcd.print("c");
         lcd.setCursor(1,0);
         lcd.print(capPointer1);
         lcd.setCursor(4,0);
         lcd.print(":");
         lcd.setCursor(5,0);
         lcd.print(capacityArray[capPointer1]);
         lcd.setCursor(0,1);
         lcd.print("c");
         lcd.setCursor(1,1);
         lcd.print(capPointer1+1);
         lcd.setCursor(4,1);
         lcd.print(":");
         lcd.setCursor(5,1);
         lcd.print(capacityArray[capPointer1+1]);
    
    if(timeElapsed3 > 5000)
     {
        timeElapsed3 = 0;
        capPointer1++;
        if(capPointer1 >= capPointer)
        {
            capPointer1=0;
        }
        
     }*/

    
    
//Send the data to ThingSpeak 
    //if(timeElapsed4 >15000)
    //{

  // ThingSpeak will only accept updates every 15 seconds.
        timeElapsed4 = 0;
    //}
    }
        lcd.setCursor(0,0);
         lcd.print("V:");
         lcd.setCursor(2,0);
         lcd.print(v);
         lcd.setCursor(0,1);
         lcd.print("I:");
         lcd.setCursor(2,1);
         lcd.print(current);
    ThingSpeak.writeField(myChannelNumber, 1, v, myWriteAPIKey);
    ThingSpeak.writeField(myChannelNumber, 2, current, myWriteAPIKey);
}
