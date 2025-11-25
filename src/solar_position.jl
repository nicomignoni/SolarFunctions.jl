using Dates, Serialization

const EARTH_PERIODIC_TERMS = open("data/earth-periodic-terms.jld", "r") do io
    deserialize(io)
end

# Limits an angle (in degrees) between 0° and 360°
mod360(angle) = mod(angle, 360)

"""
    julian_ephemeris_day(jd::Real)

A continuous count of days measured in uniform Ephemeris Time (or its successors)

# Arguments
- `jd::Real` - Julian day 
"""
function julian_ephemeris_day(jd::Real)
    ΔT = 67
    return jd + ΔT / 86400.0
end

"""
    julian_century(jd::Real)

A time interval of exactly 36,525 days (365.25 days × 100) used as a standard unit in astronomy.

# Arguments 
- `jd::Real` - Julian day
"""
function julian_century(jd::Real)
    return (jd - 2451545.0) / 36525.0 
end

"""
    julian_ephemeris_century(jed::Real)

A 36,525-day interval measured specifically in Ephemeris Time (or its modern dynamical time scales) for high-precision astronomical modeling.

# Arguments 
- `jed::Real` - Julian ephemeris day
"""
function julian_ephemeris_century(jed::Real)
    return (jed − 2451545.0) / 36525.0
end

"""
    julian_ephemeris_millenium(jec::Real)

A 1,000-year interval equal to 365,250 ephemeris days, defined within Ephemeris Time for long-term astronomical calculations.

# Arguments
- `jec::Real` - Julian ephemeris century
"""
function julian_ephemeris_millenium(jec::Real)
    return jec / 10.0
end

"""
    equation_of_time(jem::Real, Δψ::Real, ϵ::Real, α::Real)

The difference between apparent solar time and mean solar time, caused by Earth’s axial tilt and orbital eccentricity.

# Arguments 
- `jem::Real` - Julian ephemeris millenium
- `Δψ::Real` - [degrees] nutation longitude 
- `ϵ::Real` - [degrees] elliptic obliquity
- `α::Real` - [degrees] geodetic Sun's ascension 
"""
function equation_of_time(jem::Real, Δψ::Real, ϵ::Real, α::Real)
    return M - 0.0057183 - α + Δψ * cosd(ϵ)
end

function sunrise(time::DateTime)
    
end

function sunset(time::DateTime)

end

"""
    heliocentric_polynomial(Cs::Vector, jem::Real)

Polynomial approximation for the heliocentric longitude, latitude (both in radians), radius.

# Arguments
- `Cs::Vector` - Coefficients of the polynomial approximation
- `jem::Real` - Julian ephemeris millenium
"""
function heliocentric_polynomial(Cs::Vector, jem::Real)
    polynomial = 0
    for (i, C) in enumerate(Cs)
        polynomial += sum(C[:, 1] .* cos.(C[:, 2] .+ C[:, 3] .* jem)) * jem^(i - 1)
    end
    return polynomial / 1e8
end

"""
    geodetic_sun_ascension(λ::Real, ϵ::Real, β::Real)

The angle (in degrees), measured eastward from the geodetic meridian, that locates a direction or point relative to Earth’s ellipsoid rather than its rotational axis.

# Arguments
- `λ::Real` - [degrees] apparent Sun's longitude
- `ϵ::Real` - [degrees] elliptic obliquity
- `β::Real` - [degrees] geocentric latitude
"""
function geodetic_sun_ascension(λ::Real, ϵ::Real, β::Real)
    α = atand(sind(λ) * cosd(ϵ) - tand(β) * sind(ϵ), cosd(λ))
    return mod360(α)
end

"""
    geodetic_sun_declination(λ::Real, ϵ::Real, β::Real)

The Sun’s angular position north or south of Earth’s geodetic equator, measured relative to the reference ellipsoid rather than Earth’s true rotational (geocentric) equator. 

# Arguments
- `λ::Real` - [degrees] apparent Sun's longitude
- `ϵ::Real` - [degrees] elliptic obliquity
- `β::Real` - [degrees] geocentric latitude
"""
function geodetic_sun_declination(λ::Real, ϵ::Real, β::Real)
    return asind(sind(β) * cosd(ϵ) + cosd(β) * sind(ϵ) * sind(λ))
end

"""
    mean_moon_elongation_from_sun(jec::Real)

he average angular separation between the Moon and the Sun as measured along the ecliptic. 

# Arguments
- `jec::Real` - Julian ephemeris century
"""
function mean_moon_elongation_from_sun(jec::Real)
    return 297.85036 + 445267.111480jec − 0.0019142jec^2 + jec^3 / 189474.0
end

"""
    mean_sun_anomaly(jec::Real)

The angular position of the Sun in its elliptical orbit, measured from perihelion and increasing uniformly in time.

# Arguments
- `jec::Real` - Julian ephemeris century
"""
function mean_sun_anomaly(jec::Real)
    return 357.52772 + 35999.050340jec − 0.0001603jec^2 − jec^3 / 300000.0
end

"""
    mean_moon_anomaly(jec::Real)

The uniformly increasing angular position of the Moon in its elliptical orbit, measured from perigee.

# Arguments
- `jec::Real` - Julian ephemeris century
"""
function mean_moon_anomaly(jec::Real)
    return 134.96298 + 477198.867398jec + 0.0086972jec^2 + jec^3 / 56250.0
end

"""
    moon_latitude_argument(jec::Real)

The angle from the Moon’s ascending node to its position measured along its orbit, using mean (unperturbed) orbital elements.

# Arguments
- `jec::Real` - Julian ephemeris century
"""
function moon_latitude_argument(jec::Real)
    return 93.27191 + 483202.017538jec - 0.0036825jec^2 + jec^3 / 327270.0
end

"""
    moon_ascdending_longitude(jec::Real)

The ecliptic longitude, referenced to the mean equinox of the date, of the point where the Moon’s mean orbit crosses northward through the ecliptic.

# Arguments
- `jec::Real` - Julian ephemeris century
"""
function  moon_ascdending_longitude(jec::Real)
    return 125.04452 - 1934.136261jec + 0.0020708jec^2 + jec^3 / 450000.0
end

"""
    nutation_longitude(w::Vector, Ψ::Matrix, jec::Real)

The small periodic variation in Earth’s ecliptic longitude caused by gravitational torques from the Moon and Sun.

# Arguments 
- `jec::Real` - Julian ephemeris century
- `w::Vector` - [degrees] weighted nutation angles
- `Ψ::Matrix` - nutation longitude coefficients
"""
function nutation_longitude(jec::Real, w::Vector, Ψ::Matrix)
    return sum((Ψ[:, 1] .+ Ψ[:, 2] * jec) .* sind.(w)) / 3.6e7
end

"""
    nutation_obliquity(w::Vector, E::Matrix, jec::Real)

The small periodic variation in Earth’s axial tilt (obliquity) resulting from gravitational perturbations.

# Arguments
- `jec::Real` - Julian ephemeris century
- `w::Vector` - [degrees] weighted nutation angles
- `E::Matrix` - nutation obliquity coefficients
"""
function nutation_obliquity(jec::Real, w::Vector, E::Matrix)
    return sum((E[:, 1] .+ E[:, 2] * jec) .* cosd.(w)) / 3.6e7
end

"""
    elliptic_obliquity(Δϵ::Real, jem::Real)

The angle (in degrees) between Earth’s equator and the mean ecliptic defined as an ideal, unperturbed ellipse.

# Arguments: 
- `Δϵ::Real` - [degrees] nutation obliquity 
- `jem::Real` - Julian ephemeris millenium
"""
function elliptic_obliquity(jem::Real, Δϵ::Real)
    U = jem / 10.0
    ϵ₀ = 84381.448 − 4680.93U − 155.0U^2 + 1999.25U^3 − 51.38U^4 − 249.67U^5 −
         39.05U^6 + 7.12U^7 + 27.87U^8 + 5.79U^9 + 2.45U^10

    return ϵ₀ / 3600.0 + Δϵ 
end

"""
    aberration_correction(R::Real)

An adjustment applied to celestial coordinates to account for the apparent displacement caused by Earth’s motion through space.

# Arguments
- `R::Real` - [AU] heliocentric radius
"""
function aberration_correction(R::Real)
    return -20.4898 / (3600.0 * R)
end

"""
    apparent_sun_longitude(θ::Real, Δψ::Real, Δτ::Real)

The Sun’s ecliptic longitude after applying corrections for nutation and aberration.

# Arguments 
- `θ::Real` - [degrees] geocentric longitude
- `Δψ::Real` - [degrees] nutation longitude
- `Δτ::Real` - [degrees] aberration correction
"""
function apparent_sun_longitude(θ::Real, Δψ::Real, Δτ::Real)
    return θ + Δψ + Δτ
end

"""
    mean_sun_longitude(jem::Real)

The Sun’s ecliptic longitude calculated from a uniformly moving fictitious Sun on the ecliptic.

# Arguments 
- `jem::Real` - Julian ephemeris millenium
"""
function mean_sun_longitude(jem::Real)
    return 280.4664567 + 360007.6982779jem + 0.03032028jem^2 + 
           jem^3 / 49931 - jem^4 / 15300 - jem^5 / 2000000 
end

"""
    apparent_sidereal_greenwich_time(jd::Real, jc::Real, Δψ::Real, ϵ::Real)

Greenwich sidereal time corrected for the effects of nutation, giving the true rotation angle relative to the apparent equinox.

# Arguments
- `jd::Real` - Julian time
- `jc::Real` - Julian century
- `Δψ::Real` - [degress] nutation longitude 
- `ϵ::Real` - [degrees] elliptic obliquity 
"""
function apparent_sidereal_greenwich_time(jd::Real, jc::Real, Δψ::Real, ϵ::Real)
    ν₀ = mod360(280.46061837 + 360.98564736629 * (jd − 2451545.0) + 
                0.000387933jc^2 - jc^3 / 38710000.0)

    return ν₀ + Δψ * cosd(ϵ)
end

"""
    observer_local_hour_angle(ν::Real, σ::Real, α::Real)

The angle between the observer’s local meridian and a celestial object measured westward on the celestial sphere.

# Arguments
- `ν::Real` - [degrees] apparent sidereal Greenwich time 
- `σ::Real` - [degrees] observer longitude
- `α::Real` - [degrees] geodetic Sun's ascension
"""
function observer_local_hour_angle(ν::Real, σ::Real, α::Real)
    return mod360(ν + σ - α)
end

"""
    topocentric_local_hour_angle(H::Real, Δα::Real)

The hour angle of a celestial object as seen from the observer’s exact location on Earth’s surface rather than its center.

# Arguments
- `H::Real` - [degrees] observer local hour angle
- `Δα::Real` - [degrees] parallax in the Sun right ascension
"""
function topocentric_local_hour_angle(H::Real, Δα::Real)
    return H - Δα 
end

"""
    equatorial_horizontal_sun_parallax(R::Real)

The angle between the Sun’s direction as seen from Earth’s center and from a point on the equator at sea level when the Sun is on the horizon.

# Arguments 
- `R::Real` - [AU] heliocentric radius
"""
function equatorial_horizontal_sun_parallax(R::Real)
    return 8.794 / (3600.0 * R)
end

"""
    sun_ascension_parallax(x::Real, y::Real, ξ::Real, δ::Real, H::Real)

The correction to the Sun’s right ascension (in degrees) due to the difference between geocentric and topocentric perspectives.

# Arguments 
- `x::Real` - [meters] radial distance toward Earth’s equatorial plane
- `y::Real` - [meters] radial distance toward Earth’s rotational axis 
- `ξ::Real` - [degrees] equatorial horizontal parallax of the Sun
- `δ::Real` - [degrees] geodetic Sun's declination
- `H::Real` - [degrees] observer's local hour angle
"""
function sun_ascension_parallax(x::Real, y::Real, ξ::Real, δ::Real, H::Real)
    return atand(-x * sind(ξ) * sind(H), cosd(δ) - x * sind(ξ) * cosd(H))
end

"""
    topocentric_sun_ascension(α::Real, Δα::Real)

The Sun’s right ascension (in degrees) as viewed from the observer’s actual location on Earth’s surface.

# Arguments
- `α::Real` - [degrees] Sun right ascension
- `Δα::Real` - [degrees] parallax in the Sun right ascension
"""
function topocentric_sun_ascension(α::Real, Δα::Real)
    return α + Δα
end

"""
    topocentric_sun_declination(x::Real, y::Real, ξ::Real, δ::Real, H::Real, Δα::Real)

The Sun’s declination (in degrees) as viewed from the observer’s actual location on Earth’s surface.

# Arguments 
- `x::Real` - [meters] radial distance toward Earth’s equatorial plane 
- `y::Real` - [meters] radial distance toward Earth’s rotational axis
- `ξ::Real` - [degrees] equatorial horizontal parallax of the Sun
- `δ::Real` - [degrees] geodetic Sun's declination
- `H::Real` - [degrees] observer's local hour angle
- `Δα::Real` - [degrees] parallax in the Sun right ascension
"""
function topocentric_sun_declination(x::Real, y::Real, ξ::Real, δ::Real, H::Real, Δα::Real)
    return atand((sind(δ) - y * sind(ξ)) * cosd(Δα), cosd(δ) - x * sind(ξ) * cosd(H))
end

"""
    topocentric_elevation_angle(latitude::Real, temperature::Real, pressure::Real, δ̄::Real, H̄::Real)

The angle (in degrees) between the Sun and the observer’s local horizon, measured at the observer’s location.

# Arguments 
- `latitude::Real` - [degrees] observer geographical latitude
- `temperature::Real` - [Celsius] 
- `pressure::Real` - [mbars]
- `δ̄::Real` - [degrees] topocentric sun declination
- `H̄::Real` - [degrees] topocentric hour angle
"""
function topocentric_elevation_angle(
    latitude::Real,
    temperature::Real,
    pressure::Real,
    δ̄::Real,
    H̄::Real
)
    # Topocentric elevation angle without atmospheric refraction correction (in degrees)
    e₀ = asind(sind(latitude) * sind(δ̄) + cosd(latitude) * cosd(δ̄) * cosd(H̄))

    # Atmospheric refraction correction (in degrees)
    Δe =  pressure * 283.0 * 1.02 / 
          (1010.0 * (273 + temperature) * (60 * tand(e₀ + 10.3 / (e₀ + 5.11))))

    return e₀ + Δe
end

"""
    topocentric_astronomical_azimuth_angle(latitude::Real, δ̄::Real, H̄::Real)

The Sun’s azimuth (in degrees) measured from true north eastward as seen from the observer’s location, based on astronomical (not navigational) convention.

# Arguments 
- `latitude::Real` - [degrees] observer geographical latitude 
- `δ̄::Real` - [degrees] topocentric sun declination 
- `H̄::Real` - [degrees] topocentric hour angle
"""
function topocentric_astronomical_azimuth_angle(latitude::Real, δ̄::Real, H̄::Real)
    Γ = atand(sind(H̄), cosd(H̄) * sind(latitude) - tand(δ̄) * cosd(latitude)) 
    return mod360(Γ)
end

""" 
    topocentric_azimuth_angle(Γ::Real)

The Sun’s azimuth (in degrees) from the observer’s location relative to a defined reference direction—typically true north—on the horizon.

# Arguments 
- `Γ::Real` - [degrees] topocentric astronomers azimuth angle 
"""
function topocentric_azimuth_angle(Γ::Real)
    return mod360(Γ + 180) 
end

"""
    projected_incidence(e::Real, Γ::Real, ω::Real, γ::Real)

The projection (i.e., cosine) of the angle between an incoming ray (such as sunlight) and the perpendicular (normal) to the surface it strikes.

# Arguments
- `e::Real` - [degrees] Sun's elevation
- `Γ::Real` - [degrees] Sun's azimuth 
- `ω::Real` - [degrees] slope of the surface measured from the horizontal plane
- `γ::Real` - [degrees] surface azimuth rotation angle, measured from south to the projection of the surface normal on the horizontal plane, positive or negative if oriented
West or East from South, respectively.
"""
function projected_incidence(e::Real, Γ::Real, ω::Real, γ::Real)
    return sind(e) * cosd(ω) + sind(ω) * cosd(e) * cosd(Γ - γ)
end

"""
    solar_position(time::DateTime, latitude::Real, longitude::Real, altitude::Real, temperature::Real=12, pressure::Real=1013.25)

# Arguments
- `time::DateTime` - [UT] the datetime of the local observer
- `latitude::Real` - [degrees] observer geographical latitude 
- `longitude::Real` - [degrees] observer geographical longitude, positive or negative for east or west of Greenwich, respectively 
- `altitude::Real` - [meters] altitude above sea level
- `temperature::Real` - [Celsius] 
- `pressure::Real` - [mbars]

# Returns
- Sun's azimuth angle [degrees]
- Sun's elevation angle [degrees]
"""
function solar_position(
    time::DateTime, 
    latitude::Real, 
    longitude::Real,
    altitude::Real,
    temperature::Real=12,
    pressure::Real=1013.25
)

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
    θ = mod360(L + 180)

    # Geocentric latitude (in degrees)
    β = -B

    # Nutation angles (in degrees)
    X = [
        mean_moon_elongation_from_sun(jec),
        mean_sun_anomaly(jec),
        mean_moon_anomaly(jec),
        moon_latitude_argument(jec),
        moon_ascdending_longitude(jec)
    ]

    Z = EARTH_PERIODIC_TERMS.Y * X # (in degrees)

    Δψ = nutation_longitude(jec, Z, EARTH_PERIODIC_TERMS.ψ)
    Δϵ = nutation_obliquity(jec, Z, EARTH_PERIODIC_TERMS.ϵ)
    ϵ = elliptic_obliquity(jem, Δϵ)
   
    Δτ = aberration_correction(R) 
    λ = apparent_sun_longitude(θ, Δψ, Δτ)
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
