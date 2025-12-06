# API

The list above reports all the functions implemented in `SolarFunctions.jl`. In the following, we group them into topic categories and we detail them.  

```@index
```

## Solar position algorithm
The following functions implement the Solar Position Algorithm, as described in [reda2004solar](@cite). 

### Julian period conversions

```@docs
SolarFunctions.delta_T
```

```@docs
SolarFunctions.julian_ephemeris_day
```

```@docs
SolarFunctions.julian_century
```

```@docs
SolarFunctions.julian_ephemeris_century
```

```@docs
SolarFunctions.julian_ephemeris_millenium
```

### Heliocentric quantites

```@docs
SolarFunctions.heliocentric_longitude
```

```@docs
SolarFunctions.heliocentric_latitude
```

```@docs
SolarFunctions.heliocentric_radius
```

### Geocentric quantities

```@docs
SolarFunctions.geocentric_longitude
```

```@docs
SolarFunctions.geocentric_latitude
```

```@docs
SolarFunctions.geocentric_sun_ascension
```

```@docs
SolarFunctions.geocentric_sun_declination
```

### Nutation quantities

```@docs
SolarFunctions.mean_moon_elongation_from_sun
```

```@docs
SolarFunctions.mean_sun_anomaly
```

```@docs
SolarFunctions.mean_moon_anomaly
```

```@docs
SolarFunctions.moon_latitude_argument
```

```@docs
SolarFunctions.ascending_moon_longitude
```

```@docs
SolarFunctions.nutation_coefficients
```

```@docs
SolarFunctions.nutation_longitude
```

```@docs
SolarFunctions.nutation_obliquity
```

```@docs
SolarFunctions.mean_elliptic_obliquity
```

```@docs
SolarFunctions.elliptic_obliquity
```

### Solar longitude

```@docs
SolarFunctions.aberration_correction
```

```@docs
SolarFunctions.apparent_sun_longitude
```

```@docs
SolarFunctions.mean_sun_longitude
```

```@docs
SolarFunctions.mean_sidereal_greenwich_time
```

```@docs
SolarFunctions.apparent_sidereal_greenwich_time
```

```@docs
SolarFunctions.observer_local_hour
```

### Solar parallax

```@docs
SolarFunctions.reduced_observer_latitude
```

```@docs
SolarFunctions.radial_distance_equatorial_plane
```

```@docs
SolarFunctions.radial_distance_rotational_axis
```

```@docs
SolarFunctions.sun_equatorial_horizontal_parallax
```

```@docs
SolarFunctions.sun_ascension_parallax
```

### Topocentric quantities

```@docs
SolarFunctions.topocentric_sun_ascension
```

```@docs
SolarFunctions.topocentric_local_hour
```

```@docs
SolarFunctions.topocentric_sun_declination
```

```@docs
SolarFunctions.topocentric_apparent_elevation
```

```@docs
SolarFunctions.topocentric_elevation_correction
```

```@docs
SolarFunctions.topocentric_elevation
```

```@docs
SolarFunctions.topocentric_astronomical_azimuth
```

```@docs
SolarFunctions.topocentric_azimuth
```

## Clear-sky
Implements different models to calculate the GHI, DNI, and DHI components under the clear sky assumption.

```@docs
SolarFunctions.clearsky_ineichen
```

```@docs
SolarFunctions.clearsky_haurwitz
```

```@docs
SolarFunctions.clearsky_simplified_solis
```

## Atmosphere

### Linke turbidity

```@docs
SolarFunctions.aod_bb_hulstrom1980
```

```@docs
SolarFunctions.linke_turbidity_meteotest
```

```@docs
SolarFunctions.linke_turbidity_kasten1996
```

#### Altitude and pressure

```@docs
SolarFunctions.pressure2altitude
```

```@docs
SolarFunctions.altitude2pressure
```

### Airmass

```@docs
SolarFunctions.absolute_airmass
```

The following functions return the relative (not pressure-adjusted) airmass at sea level, for different models. Some models use apparent (refraction-adjusted) zenith angle while other models use true (not refraction-adjusted) zenith angle. Apparent zenith angles should be calculated at sea level. A comparison among several models is reported in [stein2012global](@cite). 

```@docs
SolarFunctions.relative_airmass_simple
```

```@docs
SolarFunctions.relative_airmass_kasten1966
```

```@docs
SolarFunctions.relative_airmass_youngirvine1967
```

```@docs
SolarFunctions.relative_arimass_kastenyoung1989
```

```@docs
SolarFunctions.relative_airmass_gueymard1993
```

```@docs
SolarFunctions.relative_airmass_young1994
```

```@docs
SolarFunctions.relative_airmass_pickering2002
```

```@docs
SolarFunctions.relative_airmass_gueymard2003
```

### Precipitable water and dew temperature

```@docs
SolarFunctions.precipitable_water_gueymard94
```

```@docs
SolarFunctions.dew_temperature2relative_humidity_oke2018
```

```@docs
SolarFunctions.relative_humidity2dew_temperature_oke2018
```

## Irradiance

```@docs
SolarFunctions.day_angle
```

```@docs
SolarFunctions.extraterrestrial_irradiance_spencer1971
```

```@docs
SolarFunctions.extraterrestrial_irradiance_asce
```

```@docs
SolarFunctions.extraterrestrial_irradiance_nrel
```
