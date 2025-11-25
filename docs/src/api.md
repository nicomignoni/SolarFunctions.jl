# API

## Irradiance
```@docs
ineichen
```

## Solar position
The following functions implement the Solar Position Algorithm, as described in [^6]. It follows the NREL's [C implementation](https://midcdmz.nrel.gov/spa/), with the only difference of using Julia's [`datetime2julian`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.julian2datetime) for the calculation of the Julian day. 

```@docs
julian_ephemeris_day
```

```@docs
julian_century
```

```@docs
julian_ephemeris_century
```

```@docs
julian_ephemeris_millenium
```

```@docs
equation_of_time
```

```@docs
heliocentric_polynomial
```

```@docs
geodetic_sun_ascension
```

```@docs
geodetic_sun_declination
```

```@docs
mean_moon_elongation_from_sun
```

```@docs
mean_sun_anomaly
```

```@docs
mean_moon_anomaly
```

```@docs
moon_latitude_argument
```

```@docs
moon_ascdending_longitude
```

```@docs
nutation_longitude
```

```@docs
nutation_obliquity
```

```@docs
elliptic_obliquity
```

```@docs
aberration_correction
```

```@docs
apparent_sun_longitude
```

```@docs
mean_sun_longitude
```

```@docs
apparent_sidereal_greenwich_time
```

```@docs
observer_local_hour_angle
```

```@docs
topocentric_local_hour_angle
```

```@docs
equatorial_horizontal_sun_parallax
```

```@docs
sun_ascension_parallax
```

```@docs
topocentric_sun_ascension
```

```@docs
topocentric_sun_declination
```

```@docs
topocentric_elevation_angle
```

```@docs
topocentric_astronomical_azimuth_angle
```

```@docs
topocentric_azimuth_angle
```

```@docs
solar_position
```

## References
[^6]: Reda, Ibrahim, and Afshin Andreas. "Solar position algorithm for solar radiation applications." Solar energy 76.5 (2004): 577-589. 
