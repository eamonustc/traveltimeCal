#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

const double pi = 3.14159265;

void directwavetime(int nl, double *velMod, double *depthMod, double delta, double depthEvent, double &takeOffAngle, double &tdir)
{
    double tolerance = 1.0e-4;
    double *thickness;
    int evtL, maxLayer, iteration;
    double distHypo, deplayerTop, thickMaxLayer;
    double lSin, uSin, uTan, lTan;
    double *aboveSinUpper, *aboveSinLower, *aboveTanUpper, *aboveTanLower;
    double deltaCalUpper, deltaCalLower, deltaDiffUpper, deltaDiffLower;
    double sinTakeOff, tanTakeOff, cosTakeOff;
    double *aboveSin, *aboveTan, *aboveCos;
    double deltaCal, deltaDiffCal;
  
    thickness = new double [nl-1];
    for (int i=0; i<nl-1; i++)
    {
        thickness[i] = depthMod[i+1] - depthMod[i];
        //std::cout << thickness[i] << std::endl;
    }

    for (int i=0; i<nl; i++)
    {
        if (depthEvent<=depthMod[i])
        {
            evtL = i-1;
            break;
        }
    }
    //std::cout << evtL << std::endl;

    if (evtL==0)
    {
        distHypo = sqrt(delta*delta + depthEvent*depthEvent);
        tdir = distHypo / velMod[0];
        takeOffAngle = 180.0 - asin(delta/distHypo)*180.0/pi;
        //std::cout << tdir << " " << takeOffAngle << std::endl;
    }
    deplayerTop = depthEvent - depthMod[evtL];

    double maxVel = 0.0;
    for (int i=0; i<=evtL; i++)
    {
        if (velMod[i] >= maxVel)
        {
            maxVel = velMod[i];
            maxLayer = i;
        }
    }
    //std::cout << maxLayer << std::endl;
    thickMaxLayer = depthMod[maxLayer+1] - depthMod[maxLayer];
    if (maxLayer==evtL)
    {
        thickMaxLayer = deplayerTop;
    }

    lSin = velMod[evtL] / maxVel * delta / sqrt(delta*delta + depthEvent*depthEvent);
    uSin = velMod[evtL] / maxVel * delta / sqrt(delta*delta + thickMaxLayer*thickMaxLayer);
    //std::cout << lSin << " " << uSin << std::endl;
    uTan = uSin / sqrt(1.0 - uSin*uSin);
    lTan = lSin / sqrt(1.0 - lSin*lSin);

    aboveSinUpper = new double [evtL];
    aboveTanUpper = new double [evtL];
    aboveSinLower = new double [evtL];
    aboveTanLower = new double [evtL];
    aboveSin = new double [evtL];
    aboveTan = new double [evtL];
    aboveCos = new double [evtL];
    deltaCalUpper = deplayerTop * uTan;
    deltaCalLower = deplayerTop * lTan;
    for (int i=0; i<evtL; i++)
    {
        aboveSinUpper[i] = velMod[i] / velMod[evtL] * uSin;
        aboveTanUpper[i] = aboveSinUpper[i] / sqrt(1.0 - aboveSinUpper[i]*aboveSinUpper[i]);
        deltaCalUpper = deltaCalUpper + thickness[i] * aboveTanUpper[i];
        aboveSinLower[i] = velMod[i] / velMod[evtL] * lSin;
        aboveTanLower[i] = aboveSinLower[i] / sqrt(1.0 - aboveSinLower[i]*aboveSinLower[i]);
        deltaCalLower = deltaCalLower + thickness[i] * aboveTanLower[i];
    }
    deltaDiffUpper = deltaCalUpper - delta;
    deltaDiffLower = deltaCalLower - delta;
    //std::cout << deltaDiffUpper << " " << deltaDiffLower << std::endl;
    
    if (std::abs(deltaCalUpper-deltaCalLower) > tolerance)
    {
        sinTakeOff = (lSin*deltaDiffUpper - uSin*deltaDiffLower) / (deltaDiffUpper - deltaDiffLower);
        tanTakeOff = sinTakeOff / sqrt(1.0 - sinTakeOff*sinTakeOff);
        deltaCal = deplayerTop * tanTakeOff;
        for (int i=0; i<evtL; i++)
        {
            aboveSin[i] = velMod[i] / velMod[evtL] * sinTakeOff;
            aboveTan[i] = aboveSin[i] / sqrt(1.0 - aboveSin[i]*aboveSin[i]);
            deltaCal = deltaCal + thickness[i] * aboveTan[i];
        }
        deltaDiffCal = deltaCal - delta;
        //std::cout << deltaDiffCal << std::endl;
        iteration = 0;

        while (std::abs(deltaDiffCal) > tolerance)
        {
            iteration = iteration + 1;
            lSin = sinTakeOff;
            lTan = lSin / sqrt(1.0 - lSin*lSin);
            deltaCalLower = deplayerTop * lTan;
            for (int i=0; i<evtL; i++)
            {
                aboveSinLower[i] = velMod[i] / velMod[evtL] * lSin;
                aboveTanLower[i] = aboveSinLower[i] / sqrt(1.0 - aboveSinLower[i]*aboveSinLower[i]);
                deltaCalLower = deltaCalLower + thickness[i] * aboveTanLower[i];
            }
            deltaDiffLower = deltaCalLower - delta;
            //std::cout << deltaDiffLower << std::endl;

            sinTakeOff = (lSin*deltaDiffUpper - uSin*deltaDiffLower) / (deltaDiffUpper - deltaDiffLower);
            tanTakeOff = sinTakeOff / sqrt(1.0 - sinTakeOff*sinTakeOff);
            deltaCal = deplayerTop * tanTakeOff;
            for (int i=0; i<evtL; i++)
            {
                aboveSin[i] = velMod[i] / velMod[evtL] * sinTakeOff;
                aboveTan[i] = aboveSin[i] / sqrt(1.0 - aboveSin[i]*aboveSin[i]);
                deltaCal = deltaCal + thickness[i] * aboveTan[i];
            }
            deltaDiffCal = deltaCal - delta;
            //std::cout << std::abs(deltaDiffCal) << " " << tolerance << std::endl;
        }
    }
    else
    {
        sinTakeOff = lSin;
        //memcpy(aboveSin, aboveSinLower, evtL*sizeof(double));
        for (int i=0; i<evtL; i++)
        {
            aboveSin[i] = aboveSinLower[i];
        }
    }
    //std::cout << sinTakeOff << std::endl;

    cosTakeOff = sqrt(1.0 - sinTakeOff*sinTakeOff);
    tdir = deplayerTop / cosTakeOff / velMod[evtL];
    for (int i=0; i<evtL; i++)
    {
        aboveCos[i] = sqrt(1.0 - aboveSin[i]*aboveSin[i]);
        tdir = tdir + thickness[i] / aboveCos[i] / velMod[i];
    }
    takeOffAngle = (pi - asin(sinTakeOff)) * 180.0 / pi;

    delete [] thickness;
    delete [] aboveSinUpper, aboveTanUpper, aboveSinLower, aboveTanLower;
    delete [] aboveSin, aboveTan, aboveCos;
}

void refractwavetime(int nl, double *velMod, double *depthMod, double delta, double depthEvent, int &refractLayer, double &takeOffAngle, double &trefract)
{
    if (nl<2)
    {
        takeOffAngle = -1;
        refractLayer = -1;
        trefract = 1000.0;
        return;
    }

    double *thickness;
    int evtL;
    double deplayerTop;
    double vRefract, maxVel;
    double sinTakeOff;

    thickness = new double [nl-1];
    for (int i=0; i<nl-1; i++)
    {
        thickness[i] = depthMod[i+1] - depthMod[i];
    }

    for (int i=0; i<nl; i++)
    {
        if (depthEvent<=depthMod[i])
        {
            evtL = i-1;
            break;
        }
    }
    deplayerTop = depthEvent - depthMod[evtL];
    //std::cout << deplayerTop << std::endl;

    double *tRefract;
    double *tid, *did;
    double *SQT, *SQTVmVn, *VnSQT;
    double *tinj, *didj;
    tRefract = new double [nl] ();
    tid = new double [nl] ();
    did = new double [nl] ();
    SQT = new double [nl] ();
    SQTVmVn = new double [nl] ();
    VnSQT = new double [nl] ();
    tinj = new double [nl] ();
    didj = new double [nl] ();
    for (int i=0; i<nl; i++)
    {
        tRefract[i] = 10000.0;
    }

    for (int m=evtL+1; m<nl; m++)
    {
        maxVel = 0.0;
        for (int i=0; i<m; i++)
        {
            if (velMod[i]>maxVel)
            {
                maxVel = velMod[i];
            }
        }
        //std::cout << maxVel << std::endl;
        if (velMod[m]>maxVel)
        {
            vRefract = velMod[m];
            did[m] = 0.0;
            tid[m] = 0.0;
            for (int i=0; i<evtL; i++)
            {
                SQT[i] = sqrt(vRefract*vRefract - velMod[i]*velMod[i]);
                SQTVmVn[i] = SQT[i] / vRefract / velMod[i];
                VnSQT[i] = velMod[i] / SQT[i];
                did[m] = did[m] + thickness[i] * VnSQT[i];
                tid[m] = tid[m] + thickness[i] * SQTVmVn[i];
            }
            for (int i=evtL; i<m; i++)
            {
                SQT[i] = sqrt(vRefract*vRefract - velMod[i]*velMod[i]);
                SQTVmVn[i] = SQT[i] / vRefract / velMod[i];
                VnSQT[i] = velMod[i] / SQT[i];
                did[m] = did[m] + 2.0 * thickness[i] * VnSQT[i];
                tid[m] = tid[m] + 2.0 * thickness[i] * SQTVmVn[i];
            }
            didj[m] = did[m] - deplayerTop*VnSQT[evtL];
            tinj[m] = tid[m] - deplayerTop*SQTVmVn[evtL];

            if (didj[m]<delta)
            {
                tRefract[m] = tinj[m] + delta / vRefract;
            }
            else
            {
                tRefract[m] = 10000.0;
            }
        }
        else
        {
            tid[m] = 10000.0;
            did[m] = 10000.0;
            tRefract[m] = 10000.0;
        }
    }

    trefract = 10000.0;
    for (int i=0; i<nl; i++)
    {
        if (trefract>tRefract[i])
        {
            trefract = tRefract[i];
            refractLayer = i;
        }
    }
    if (trefract<10000.0)
    {
        sinTakeOff  = velMod[evtL] / velMod[refractLayer];
        takeOffAngle = asin(sinTakeOff) * 180.0 / pi;
    }
    else
    {
        takeOffAngle = -1;
        refractLayer = -1;
        trefract = 10000.0;
    }

    delete [] thickness;
    delete [] tRefract;
    delete [] tid, did;
    delete [] SQT, SQTVmVn, VnSQT;
    delete [] tinj, didj;
}

void traveltimeCal(int nl, double *velMod, double *depthMod, double delta, double depthEvent, double &takeOffAngle, double &travelTime, int &isRefract, int &refractLayer)
{
    double takeOffAngleDir, takeOffAngleRef;
    double tdir, trefract;

    isRefract = 0;
    refractwavetime(nl, velMod, depthMod, delta, depthEvent, refractLayer, takeOffAngleRef, trefract);
    directwavetime(nl, velMod, depthMod, delta, depthEvent, takeOffAngleDir, tdir);
    if (tdir<trefract)
    {
        travelTime = tdir;
        takeOffAngle = takeOffAngleDir;
        refractLayer = -1;
    }
    else
    {
        travelTime = trefract;
        takeOffAngle = takeOffAngleRef;
        isRefract = 1;
    }
}

int main()
{
    int nl = 9;
    int isRefract, refractLayer;
    double takeOffAngle, travelTime;
    double velMod[9] = {4.5, 5.5, 6.75, 6.75, 6.75, 6.9, 7.75, 8.0, 8.175};
    double depthMod[9] = {0.0, 4.0, 10.0, 20.0, 30.0, 35.0, 40.0, 150.0, 165.0};
    double delta = 100.0;
    double depthEvent = 6.5;

    traveltimeCal(nl, velMod, depthMod, delta, depthEvent, takeOffAngle, travelTime, isRefract, refractLayer);
    std::cout << " takeOffAngle: " << takeOffAngle << " travelTime: " << travelTime << std::endl;
}
