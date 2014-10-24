/*
 * Class to represent a geographic position given by a latitude and longitude.
 *
 * 2013-05-28   Geoff Hayes     Initial Release.
 */
#include "Position.h"
#include <cmath>

Position::Position(double lat, double lon, VALUE_TYPE units)
{
    switch (units)
    {
        case RADIANS:
        {
            _latitude  = lat;
            _longitude = lon;
            break;
        }
        case DEGREES:
        {
            _latitude  = lat*M_PI/180.0;
            _longitude = lon*M_PI/180.0;
            break;
        }
    }

    // need to ensure that the latitude and longitude are bounded
    if (_latitude > M_PI/2.0)
    {
        _latitude = M_PI/2.0;
    }
    else if (_latitude < -M_PI/2.0)
    {
        _latitude = -M_PI/2.0;
    }

    if (_longitude > M_PI)
    {
        _longitude = M_PI;
    }
    else if (_longitude < -M_PI)
    {
        _longitude = -M_PI;
    }
}

Position::~Position()
{
    // intentionally left blank
}

double Position::getLatitude(VALUE_TYPE unit) const
{
    double lat = _latitude;

    switch (unit)
    {
        case DEGREES:
        {
            lat = _latitude*180.0/M_PI;
            break;
        }
        case RADIANS:
        {
            // nothing to do
            break;
        }
    }

    return lat;
}

double Position::getLongitude(VALUE_TYPE unit) const
{
    double lon = _longitude;

    switch (unit)
    {
        case DEGREES:
        {
            lon = _longitude*180.0/M_PI;
            break;
        }
        case RADIANS:
        {
            // nothing to do
            break;
        }
    }

    return lon;
}

void Position::setLatitude(double lat, VALUE_TYPE unit)
{
    switch (unit)
    {
        case DEGREES:
        {
            _latitude = lat*180.0/M_PI;
            break;
        }
        case RADIANS:
        {
            _latitude = lat;
            break;
        }
    }

    // ensure that the latitude is in the correct range
    if (_latitude > M_PI/2.0)
    {
        _latitude = M_PI/2.0;
    }
    else if (_latitude < -M_PI/2.0)
    {
        _latitude = -M_PI/2.0;
    }
}

void Position::setLongitude(double lon, VALUE_TYPE unit)
{
    switch (unit)
    {
        case DEGREES:
        {
            _longitude = lon*180.0/M_PI;
            break;
        }
        case RADIANS:
        {
            _longitude = lon;
            break;
        }
    }

    // ensure that the latitude is in the correct range
    if (_longitude > M_PI)
    {
        _longitude = M_PI;
    }
    else if (_longitude < -M_PI)
    {
        _longitude = -M_PI;
    }
}




