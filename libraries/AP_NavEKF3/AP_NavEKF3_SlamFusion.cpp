/// -*- tab-width: 4; Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*-

#include <AP_HAL/AP_HAL.h>

#if HAL_CPU_CLASS >= HAL_CPU_CLASS_150

#include "AP_NavEKF3.h"
#include "AP_NavEKF3_core.h"
#include <AP_AHRS/AP_AHRS.h>
#include <AP_Vehicle/AP_Vehicle.h>

extern const AP_HAL::HAL& hal;


/********************************************************
*                   FUSE MEASURED_DATA                  *
********************************************************/
// select fusion of position estimates from SLAM system
void NavEKF3_core::SelectSlamFusion()
{
    // Check if the magnetometer has been fused on that time step and the filter is running at faster than 200 Hz
    // If so, don't fuse measurements on this time step to reduce frame over-runs
    // Only allow one time slip to prevent high rate magnetometer data preventing fusion of other measurements
    if (magFusePerformed && dtIMUavg < 0.005f && !posVelFusionDelayed) {
        slamFusionDelayed = true;
        return;
    } else {
        slamFusionDelayed = false;
    }

    // read GPS data from the sensor and check for new data in the buffer
    readSlamData();
    slamDataToFuse = storedSLAM.recall(slamDataDelayed,imuDataDelayed.time_ms);
    
    // perform fusion
    if (slamDataToFuse) {
        FuseSlamNED(); 
    }
}

// fuse selected position and height measurements
void NavEKF3_core::FuseSlamNED()
{    
    Vector3f Innov;          
    Vector3f observation;    


    // form the observation vector
    observation[0] = slamDataDelayed.pos.x;
    observation[1] = slamDataDelayed.pos.y;
    observation[2] = slamDataDelayed.pos.z;

    // Compute innovation
    Innov[0] = stateStruct.position.x - observation[0];
    Innov[1] = stateStruct.position.y - observation[1];
    Innov[2] = stateStruct.position.z - observation[2];
    
    // Observation covariance
    Matrix3f R_OBS(slamDataDelayed.posCov[0],slamDataDelayed.posCov[1],slamDataDelayed.posCov[2],
                   slamDataDelayed.posCov[3],slamDataDelayed.posCov[4],slamDataDelayed.posCov[5],
                   slamDataDelayed.posCov[6],slamDataDelayed.posCov[7],slamDataDelayed.posCov[8]);

    // Covariance of mean position estimate
    Matrix3f Sigma(P[7][7], P[7][8], P[7][9],
                   P[8][7], P[8][8], P[8][9],
                   P[9][7], P[9][8], P[9][9]);
    
    // Note: Measurement jacobian is Identity just like the GPS fusion.
    Matrix3f varInnov = Sigma + R_OBS;
    Matrix3f invVarInnov;
    bool inverseExists = varInnov.inverse(invVarInnov); 

    if(!inverseExists){    
        return;
    }

    Vector24 sKfusion1,sKfusion2,sKfusion3; // Each column of the Kalman gain matrix
    for(uint8_t i = 0; i<= stateIndexLim; i++){            
        sKfusion1[i] = P[i][7]*invVarInnov[0][0] + P[i][8]*invVarInnov[1][0] + P[i][9]*invVarInnov[2][0];
        sKfusion2[i] = P[i][7]*invVarInnov[0][1] + P[i][8]*invVarInnov[1][1] + P[i][9]*invVarInnov[2][1];
        sKfusion3[i] = P[i][7]*invVarInnov[0][2] + P[i][8]*invVarInnov[1][2] + P[i][9]*invVarInnov[2][2];
    }

    // Efficient matrix multiplication of K*H*P
    Matrix24 sKHP;
    for (uint8_t i= 0; i<=stateIndexLim; i++) {
        for (uint8_t j= 0; j<=stateIndexLim; j++)
            {
               sKHP[i][j] = sKfusion1[i] * P[7][j];
            }
    }

    for (uint8_t i= 0; i<=stateIndexLim; i++) {
        for (uint8_t j= 0; j<=stateIndexLim; j++)
            {
               sKHP[i][j] += sKfusion2[i] * P[8][j];
            }
    }

    for (uint8_t i= 0; i<=stateIndexLim; i++) {
        for (uint8_t j= 0; j<=stateIndexLim; j++)
            {
               sKHP[i][j] += sKfusion3[i] * P[9][j];
            }
    }
        
    // Check that we are not going to drive any variances negative and skip the update if so
    bool healthyFusion = true;
    for (uint8_t i= 0; i<=stateIndexLim; i++) {
        if (sKHP[i][i] > P[i][i]) {
            healthyFusion = false;
        }
    }
    
    if(healthyFusion){
        // update the covariance matrix
        for (uint8_t i= 0; i<=stateIndexLim; i++) {
            for (uint8_t j= 0; j<=stateIndexLim; j++) {
                P[i][j] = P[i][j] - KHP[i][j];
            }
        }

        // force the covariance matrix to be symmetrical and limit the variances to prevent ill-condiioning.
        ForceSymmetry();
        ConstrainVariances();
        
        // update states and renormalise the quaternions
        for (uint8_t i = 0; i<=stateIndexLim; i++) {
            statesArray[i] = statesArray[i] - sKfusion1[i] * Innov[0] - sKfusion2[i] *Innov[1] - sKfusion3[i]*Innov[2] ;
        }

        stateStruct.quat.normalize();                
    }            
}



#endif // HAL_CPU_CLASS
