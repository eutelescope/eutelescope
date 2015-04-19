#ifndef EUTELGEOMETRYTELECOPEGEODESCRIPTION_TCC  
#define EUTELGEOMETRYTELECOPEGEODESCRIPTION_TCC

/** Determine id of the sensor in which point is locate
 *  * 
 *  * @param globalPos 3D point in global reference frame
 *  * @return sensorID or -999 if the point in outside of sensor volume
 *  */
//MUST OUTPUT -999 TO SIGNIFY THAT NO SENSOR HAS BEEN FOUND. SINCE USED IN PATTERN RECOGNITION THIS WAY.
template<class number>
int EUTelGeometryTelescopeGeoDescription::getSensorID( number globalPos[] ) const {
    streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::getSensorID() " << std::endl;
    const float constPos[3] = {globalPos[0],globalPos[1],globalPos[2]};
    _geoManager->FindNode( constPos[0], constPos[1], constPos[2] );

    std::vector<std::string> split;

    int sensorID = -999;

    const char* volName1 = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
    streamlog_out(DEBUG2) << "init sensorID  : " << sensorID  <<  " " << volName1 << std::endl;

    while( _geoManager->GetLevel() > 0 ) {
        const char* volName = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
        streamlog_out( DEBUG1 ) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") found in volume: " << volName << " level: " << _geoManager->GetLevel() << std::endl;
        split = Utility::stringSplit( std::string( volName ), "/", false);
        if ( split.size() > 0 && split[0].length() > 16 && (split[0].substr(0,16) == "volume_SensorID:") ) {
            int strLength = split[0].length();
            sensorID = strtol( (split[0].substr(16, strLength )).c_str(), NULL, 10 );
            streamlog_out(DEBUG1) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was found at :" << sensorID << std::endl;
        break;
        }
    _geoManager->CdUp();  ////////////////////////////////////////THIS NEEDS TO BE FIXED. If partice falls in the pixel volume and to find sensor ID you need to be on the sensor volume
    }

    const char* volName2 = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
    streamlog_out( DEBUG2 ) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") found in volume: " << volName2 << " no moving around any more" << std::endl;

    if( sensorID >= 0 )
    {
        //                sensorID = strtol( (split[0].substr(16, strLength )).c_str(), NULL, 10 );
        streamlog_out(DEBUG5) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was found. sensorID = " << sensorID << std::endl;
    }
    else
    {
        streamlog_out(DEBUG5) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was not found inside any sensor! sensorID = " << sensorID << std::endl;
    }
    return sensorID;
}
#endif
