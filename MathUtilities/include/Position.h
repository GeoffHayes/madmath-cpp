/*
 * Class to represent a geographic position given by a latitude and longitude.
 *
 * 2013-05-28   Geoff Hayes     Initial Release.
 */

#ifndef POSITION_H_
#define POSITION_H_

//! Class to represent a geographic position given by a latitude and longitude.
class Position
{
    public:

        //! Enum to indicate the inputs or output units to an instance
        //! of this class.
        enum VALUE_TYPE
        {
            RADIANS=0, /**< enumeration for input/output as radians */
            DEGREES=1  /**< enumeration for input/output as degrees */
        };

        /**
         * Class constructor to initialize the latitude and longitude members.
         *
         * @param   lat   Position latitude in units matching input type
         *                (default value is zero).
         * @param   lon   Position longitude in units matching input type
         *                (default value is zero).
         * @param   unit  Units of input (default is radians).
         */
        Position(double lat=0.0, double lon=0.0, VALUE_TYPE unit=RADIANS);

        /**
         * Class destructor
         */
        ~Position();

        /**
         * Returns the latitude in the requested units.
         *
         * @param   unit   Units of requested latitude (default is RADIANS).
         *
         * @return  The latitude in the requested units.
         */
        double getLatitude(VALUE_TYPE unit=RADIANS) const;

        /**
         * Returns the longitude in the requested units.
         *
         * @param   unit   Units of requested longitude (default is RADIANS).
         *
         * @return  The longitude in the requested units.
         */
        double getLongitude(VALUE_TYPE unit=RADIANS) const;

        /**
         * Sets the latitude in the requested units.
         *
         * @param   lat    Postion latitude to be set.
         * @param   unit   Units of requested latitude (default is RADIANS).
         */
        void setLatitude(double lat, VALUE_TYPE unit=RADIANS);

        /**
         * Returns the longitude in the requested units.
         *
         * @param   lon    Position longitude to be set.
         * @param   unit   Units of requested longitude (default is RADIANS).
         *
         */
        void setLongitude(double lon, VALUE_TYPE unit=RADIANS);

    private:

        //! Position latitude in radians (native type) bounded by [-pi/2, pi/2]
        double _latitude;

        //! Position longitude in radians (native type) bounded by [-pi, pi]
        double _longitude;
};




#endif /* POSITION_H_ */
