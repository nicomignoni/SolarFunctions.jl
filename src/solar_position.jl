using Dates, Serialization

const EARTH_PERIODIC_TERMS = open("data/earth-periodic-terms.jld", "r") do io
    deserialize(io)
end

# Limits an angle (in degrees) between 0° and 360°
mod360(angle) = mod(angle, 360)

function julian_ephemeris_day(jd::Real)
    """
        julian_ephemeris_day(jd::Real)

    A continuous count of days measured in uniform Ephemeris Time (or its successors)

    # Arguments
    - `jd::Real` - Julian day 
    """
    ΔT = 67
    return jd + ΔT / 86400.0
end

function julian_century(jd::Real)
    """
        julian_century(jd::Real)

    A time interval of exactly 36,525 days (365.25 days × 100) used as a standard unit in astronomy.

    # Arguments 
    - `jd::Real` - Julian day
    """
    return (jd - 2451545.0) / 36525.0 
end

function julian_ephemeris_century(jed::Real)
    """
        julian_ephemeris_century(jed::Real)

    A 36,525-day interval measured specifically in Ephemeris Time (or its modern dynamical time scales) for high-precision astronomical modeling.

    # Arguments 
    - `jed::Real` - Julian ephemeris day
    """
    return (jed − 2451545.0) / 36525.0
end

function julian_ephemeris_millenium(jec::Real)
    """
        julian_ephemeris_millenium(jec::Real)

    A 1,000-year interval equal to 365,250 ephemeris days, defined within Ephemeris Time for long-term astronomical calculations.

    # Arguments
    - `jec::Real` - Julian ephemeris century
    """
    return jec / 10.0
end

function heliocentric_polynomial(Cs::Vector, jem::Real)
    """
        heliocentric_polynomial(Cs::Vector, jem::Real)

    Polynomial approximation for the heliocentric longitude, latitude (both in radians), radius.

    # Arguments
    - `Cs::Vector` - Coefficients of the polynomial approximation
    - `jem::Real` - Julian ephemeris millenium
    """
    polynomial = 0
    for (i, C) in enumerate(Cs)
        polynomial += sum(C[:, 1] .* cos.(C[:, 2] .+ C[:, 3] .* jem)) * jem^(i - 1)
    end
    return 1e-8polynomial 
end

function geodetic_sun_ascension(λ::Real, ϵ::Real, β::Real)
    """
        geodetic_sun_ascension(λ::Real, ϵ::Real, β::Real)

    The angle (in degrees), measured eastward from the geodetic meridian, that locates a direction or point relative to Earth’s ellipsoid rather than its rotational axis.

    # Arguments
    - `λ::Real` - [degrees] apparent Sun's longitude
    - `ϵ::Real` - [degrees] elliptic obliquity
    - `β::Real` - [degrees] geocentric latitude
    """
    return atand(sind(λ) * cosd(ϵ) - tand(β) * sind(ϵ), cosd(λ))
end

function geodetic_sun_declination(λ::Real, ϵ::Real, β::Real)
    """
        geodetic_sun_declination(λ::Real, ϵ::Real, β::Real)

    The Sun’s angular position north or south of Earth’s geodetic equator, measured relative to the reference ellipsoid rather than Earth’s true rotational (geocentric) equator. 

    # Arguments
    - `λ::Real` - [degrees] apparent Sun's longitude
    - `ϵ::Real` - [degrees] elliptic obliquity
    - `β::Real` - [degrees] geocentric latitude
    """

    return asind(sind(β) * cosd(ϵ) + cosd(β) * sind(ϵ) * sind(λ))
end

function mean_moon_elongation_from_sun(jec::Real)
    """
        mean_moon_elongation_from_sun(jec::Real)

   The average angular separation between the Moon and the Sun as measured along the ecliptic. 

    # Arguments
    - `jec::Real` - Julian ephemeris century
    """
    return 297.85036 + 445267.111480 * jec − 0.0019142 * jec^2 + jec^3 / 189474.0
end

function mean_sun_anomaly(jec::Real)
    """
        mean_sun_anomaly(jec::Real)

    The angular position of the Sun in its elliptical orbit, measured from perihelion and increasing uniformly in time.

    # Arguments
    - `jec::Real` - Julian ephemeris century
    """
    return 357.52772 + 35999.050340 * jec − 0.0001603 * jec^2 − jec^3 / 300000.0
end

function mean_moon_anomaly(jec::Real)
    """
        mean_moon_anomaly(jec::Real)

    The uniformly increasing angular position of the Moon in its elliptical orbit, measured from perigee.

    # Arguments
    - `jec::Real` - Julian ephemeris century
    """
    return 134.96298 + 477198.867398 * jec + 0.0086972 * jec^2 + jec^3 / 56250.0
end

function moon_latitude_argument(jec::Real)
    """
        moon_latitude_argument(jec::Real)

    The angle from the Moon’s ascending node to its position measured along its orbit, using mean (unperturbed) orbital elements.

    # Arguments
    - `jec::Real` - Julian ephemeris century
    """
    return 93.27191 + 483202.017538 * jec - 0.0036825 * jec^2 + jec^3 / 327270.0
end

function  moon_ascdending_longitude(jec::Real)
    """
        moon_ascdending_longitude(jec::Real)

    The ecliptic longitude, referenced to the mean equinox of the date, of the point where the Moon’s mean orbit crosses northward through the ecliptic.

    # Arguments
    - `jec::Real` - Julian ephemeris century
    """
    return 125.04452 - 1934.136261 * jec + 0.0020708 * jec^2 + jec^3 / 450000.0
end

function nutation_longitude(jec::Real, w::Vector, Ψ::Matrix)
    """
        nutation_longitude(w::Vector, Ψ::Matrix, jec::Real)

    The small periodic variation in Earth’s ecliptic longitude caused by gravitational torques from the Moon and Sun.

    # Arguments 
    - `w::Vector` - [degrees] weighted nutation angles
    - `Ψ::Matrix` - nutation longitude coefficients
    - `jec::Real` - Julian ephemeris century
    """
    return 3.6e-7sum((Ψ[:, 1] .+ Ψ[:, 2] * jec) .* sind.(w))
end

function nutation_obliquity(jec::Real, w::Vector, E::Matrix)
    """
        nutation_obliquity(w::Vector, E::Matrix, jec::Real)

    The small periodic variation in Earth’s axial tilt (obliquity) resulting from gravitational perturbations.

    # Arguments
    - `w::Vector` - [degrees] weighted nutation angles
    - `E::Matrix` - nutation obliquity coefficients
    - `jec::Real` - Julian ephemeris century
    """
    return 3.6e-7sum((E[:, 1] .+ E[:, 2] * jec) .* cosd.(w))
end

function elliptic_obliquity(jem::Real, Δϵ::Real)
    """
        elliptic_obliquity(Δϵ::Real, jem::Real)

    The angle (in degrees) between Earth’s equator and the mean ecliptic defined as an ideal, unperturbed ellipse.

    # Arguments: 
    - `Δϵ::Real` - [degrees] nutation obliquity 
    - `jem::Real` - Julian ephemeris millenium
    """
    U = jem / 10.0
    ϵ₀ = 84381.448 − 4680.93U − 155.0U^2 + 1999.25U^3 − 51.38U^4 − 249.67U^5 −
         39.05U^6 + 7.12U^7 + 27.87U^8 + 5.79U^9 + 2.45U^10

    return ϵ₀ / 3600.0 + Δϵ 
end

function aberration_correction(R::Real)
    """
        aberration_correction(R::Real)

    # Arguments
    - `R::Real` - [AU] heliocentric radius
    """
    return -20.4898 / (3600.0 * R)
end

function apparent_sun_longitude(Θ::Real, Δψ::Real, Δτ::Real)
    """
        apparent_sun_longitude(Θ::Real, Δψ::Real, Δτ::Real)

    # Arguments 
    - `Θ::Real` - [degrees] geocentric longitude
    - `Δψ::Real` - [degrees] nutation longitude
    - `Δτ::Real` - [degrees] aberration correction
    """
    return Θ + Δψ + Δτ
end

function mean_sun_longitude(jem::Real)
    """
        mean_sun_longitude(jem::Real)

    # Arguments 
    - `jem::Real` - Julian ephemeris millenium
    """
    return 280.4664567 + 360007.6982779jem + 0.03032028jem^2 + 
           jem^3 / 49931 - jem^4 / 15300 - jem^5 / 2000000 
end

function apparent_sidereal_greenwich_time(jd::Real, jc::Real, Δψ::Real, ϵ::Real)
    """
        apparent_sidereal_greenwich_time(jd::Real, jc::Real, Δψ::Real, ϵ::Real)

    # Arguments
    - `jd::Real` - Julian time
    - `jc::Real` - Julian century
    - `Δψ::Real` - nutation longitude 
    - `ϵ::Real` - elliptic obliquity 
    """
    ν₀ = mod360(280.46061837 + 360.98564736629 * (jd − 2451545.0) + 
                0.000387933jc^2 - jc^3 / 38710000.0)

    return ν₀ + Δψ * cosd(ϵ)
end

function observer_local_hour_angle(ν::Real, σ::Real, α::Real)
    """
        observer_local_hour_angle(ν::Real, σ::Real, α::Real)

    # Arguments
    - `ν::Real` - [degrees] apparent sidereal Greenwich time 
    - `σ::Real` - [degrees] observer longitude
    - `α::Real` - [degrees] geodetic Sun's ascension
    """
    return mod360(ν + σ - α)
end

function topocentric_local_hour_angle(H::Real, Δα::Real)
    """
        topocentric_local_hour_angle(H::Real, Δα::Real)

    # Arguments
    - `H::Real` - [degrees] observer local hour angle
    - `Δα::Real` - [degrees] parallax in the Sun right ascension
    """
    return H - Δα 
end

function equatorial_horizontal_sun_parallax(R::Real)
    """
        equatorial_horizontal_sun_parallax(R::Real)

    # Arguments 
    - `R::Real` - [AU] heliocentric radius
    """
    return 8.794 / (3600.0 * R)
end

function sun_ascension_parallax(x::Real, y::Real, ξ::Real, δ::Real, H::Real)
    """
        sun_ascension_parallax(x::Real, y::Real, ξ::Real, δ::Real, H::Real)

    # Arguments 
    - `x::Real` - [meters] radial distance toward Earth’s equatorial plane
    - `y::Real` - [meters] radial distance toward Earth’s rotational axis 
    - `ξ::Real` - [degrees] equatorial horizontal parallax of the Sun
    - `δ::Real` - [degrees] geodetic Sun's declination
    - `H::Real` - [degrees] observer's local hour angle
    """
    return atand(-x * sind(ξ) * sind(H), cosd(δ) - x * sind(ξ) * cosd(H))
end

function topocentric_sun_ascension(α::Real, Δα::Real)
    """
        topocentric_sun_ascension(α::Real, Δα::Real)

    # Arguments
    - `α::Real` - [degrees] Sun right ascension
    - `Δα::Real` - [degrees] parallax in the Sun right ascension
    """
    return α + Δα
end

function topocentric_sun_declination(x::Real, y::Real, ξ::Real, δ::Real, H::Real, Δα::Real)
    """
        topocentric_sun_declination(x::Real, y::Real, ξ::Real, δ::Real, H::Real, Δα::Real)

    # Arguments 
    - `x::Real` - [meters] radial distance toward Earth’s equatorial plane 
    - `y::Real` - [meters] radial distance toward Earth’s rotational axis
    - `ξ::Real` - [degrees] equatorial horizontal parallax of the Sun
    - `δ::Real` - [degrees] geodetic Sun's declination
    - `H::Real` - [degrees] observer's local hour angle
    - `Δα::Real` - [degrees] parallax in the Sun right ascension
    """
    return atand((sind(δ) - y * sind(ξ)) * cosd(Δα), cosd(δ) - x * sind(ξ) * cosd(H))
end

function topocentric_elevation_angle(
    latitude::Real,
    temperature::Real,
    pressure::Real,
    δ̄::Real,
    H̄::Real
)
    """
        topocentric_elevation_angle(latitude::Real, temperature::Real, pressure::Real, δ̄::Real, H̄::Real)

    # Arguments 
    - `latitude::Real` - [degrees] observer geographical latitude
    - `temperature::Real` - [Celsius] 
    - `pressure::Real` - [mbars]
    - `δ̄::Real` - [degrees] topocentric sun declination
    - `H̄::Real` - [degrees] topocentric hour angle
    """
    # Topocentric elevation angle without atmospheric refraction correction (in degrees)
    e₀ = asind(sind(latitude) * sind(δ̄) + cosd(latitude) * cosd(δ̄) * cosd(H̄))

    # Atmospheric refraction correction (in degrees)
    Δe =  pressure * 283.0 * 1.02 / 
          (1010.0 * (273 + temperature) * (60 * tand(e₀ + 10.3 / (e₀ + 5.11))))

    return e₀ + Δe

end

function topocentric_astronomical_azimuth_angle(latitude::Real, δ̄::Real, H̄::Real)
    """
        topocentric_astronomical_azimuth_angle(latitude::Real, δ̄::Real, H̄::Real)

    # Arguments 
    - `latitude::Real` - [degrees] observer geographical latitude 
    - `δ̄::Real` - [degrees] topocentric sun declination 
    - `H̄::Real` - [degrees] topocentric hour angle
    """
    Γ = atand(sind(H̄), cosd(H̄) * sind(latitude) - tand(δ̄) * cosd(latitude)) 
    return mod360(Γ)
end

function topocentric_azimuth_angle(Γ::Real)
    """ 
        topocentric_azimuth_angle(Γ::Real)

    # Arguments 
    - `Γ::Real` - [degrees] topocentric astronomers azimuth angle 
    """
    return mod360(Γ + 180) 
end

function solar_position(
    time::DateTime, 
    latitude::Real, 
    longitude::Real,
    altitude::Real,
    temperature::Real=12,
    pressure::Real=1013.25
)
    """
        solar_position(time::DateTime, latitude::Real, longitude::Real, altitude::Real, temperature::Real=12, pressure::Real=1013.25)

    # Arguments
    - `time::DateTime` - [UT] the datetime of the local observer
    - `latitude::Real` - [degrees] observer geographical latitude 
    - `longitude::Real` - [degrees] observer geographical longitude, positive or negative for east or west of Greenwich, respectively 
    - `altitude::Real` - [meters]
    - `temperature::Real` - [Celsius]
    - `pressure::Real` - [mbars]

    # Returns
    - Sun's azimuth angle [degrees]
    - Sun's elevation angle [degrees]
    """

    jd = datetime2julian(time)
    jed = julian_ephemeris_day(jd) 
    jc = julian_century(jd) 
    jec = julian_ephemeris_century(jed)
    jem = julian_ephemeris_millenium(jec) 

    # Heliocentric longitude (in degrees)
    L = heliocentric_polynomial(EARTH_PERIODIC_TERMS.L, jem) |> rad2deg |> mod360

    # Heliocentric latitude (in degrees) 
    B = heliocentric_polynomial(EARTH_PERIODIC_TERMS.B, jem) |> rad2deg
    
    # Heliocentric radius (in astronomical units [AU])
    R = heliocentric_polynomial(EARTH_PERIODIC_TERMS.R, jem)

    # Geocentric longitude (in degrees)
    Θ = mod360(L + 180)

    # Geocentric latitude (in degrees)
    β = -B

    # Nutation angles (in degrees)
    X = [
        mean_moon_elongation_from_sun(jec)
        mean_sun_anomaly(jec)
        mean_moon_anomaly(jec)
        moon_latitude_argument(jec)
        moon_ascdending_longitude(jec)
    ]

    Z = EARTH_PERIODIC_TERMS.Y * X # (in degrees)

    Δψ = nutation_longitude(jec, Z, EARTH_PERIODIC_TERMS.ψ)
    Δϵ = nutation_obliquity(jec, Z, EARTH_PERIODIC_TERMS.ψ)
    ϵ = elliptic_obliquity(jem, Δϵ)
   
    Δτ = aberration_correction(R) 
    λ = apparent_sun_longitude(Θ, Δψ, Δτ)
    ν = apparent_sidereal_greenwich_time(jd, jc, Δψ, ϵ) 

    α = geodetic_sun_ascension(λ, ϵ, β)
    δ = geodetic_sun_declination(λ, ϵ, β)
    H = observer_local_hour_angle(ν, longitude, α)
    ξ = equatorial_horizontal_sun_parallax(R)
    
    u = atand(0.99664719 * tand(latitude)) # (in degrees)
    x = cosd(u) + altitude * cosd(latitude) / 6378140.0 
    y = 0.99664719 * sind(u) + altitude * sind(latitude) / 6378140.0

    Δα = sun_ascension_parallax(x, y, ξ, δ, H) 

    δ̄ = topocentric_sun_declination(x, y, ξ, δ, H, Δα)
    ᾱ = topocentric_sun_ascension(α, Δα)
    H̄ = topocentric_local_hour_angle(H, Δα)
    
    e = topocentric_elevation_angle(latitude, temperature, pressure, δ̄, H̄)

    Γ = topocentric_astronomical_azimuth_angle(latitude, δ̄, H̄)
    Φ = topocentric_azimuth_angle(Γ)

    return Φ, e
end

function equation_of_time(jem::Real, Δψ::Real, ϵ::Real, α::Real)
    """
    # Arguments 
    - `jem::Real` - Julian ephemeris millenium
    - `Δψ::Real` - 
    - `ϵ::Real` - 
    - `α::Real` - 
    """
    return M - 0.0057183 - α + Δψ * cosd(ϵ)
end

function sunrise(time::DateTime)
    
end

function sunset(time::DateTime)

end
