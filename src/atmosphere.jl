using Pkg.Artifacts, Serialization

# The linke-turbidity.jld contains a single 2160 x 4320 x 12 array of type uint8. 
# To determine the Linke urbidity for a position on the Earth's surface 
# for a given month do the following: LT = LinkeTurbidity(LatitudeIndex, LongitudeIndex, month).

# The nodes of the grid are 5' (1/12=0.0833[arcdeg]) apart.
# From Section 8 of Aerosol optical depth and Linke turbidity climatology
# http://www.meteonorm.com/images/uploads/downloads/ieashc36_report_TL_AOD_climatologies.pdf
# 1st row: 89.9583 S, 2nd row: 89.875 S
# 1st column: 179.9583 m^2, 2nd column: 179.875 m^2
const LINKE_TURBIDITY_METEOTEST = open(artifact"data/linke-turbidity.jld", "r") do io;
    deserialize(io)
end

# Atmospheric pressure at sea level [Pa]
const ATMOSPHERIC_PRESSURE = 101325.0

# Bounds [deg] for latitude and longitude
const LATITUDE_MAX = 90
const LATITUDE_MIN = -90
const LONGITUDE_MAX = 180
const LONGITUDE_MIN = -180 

const LATITUDE_RANGE = LATITUDE_MAX - LATITUDE_MIN
const LONGITUDE_RANGE = LONGITUDE_MAX - LONGITUDE_MIN

"""
    aod_bb_hulstrom1980(aod380::Real, aod500::Real)

Approximate broadband aerosol optical depth (AOD) as proposed in [bird1980direct, bird1981review](@cite), which develops a correlation for broadband aerosol optical depth using two wavelengths, 380 nm and 500 nm. Arguments `aod380` and `aod500` are the AOD measured at 380 nm and 500 nm, respectively.
"""
function aod_bb_hulstrom1980(aod380::Real, aod500::Real)
    return 0.27583aod380 + 0.35aod500
end

"""
    linke_turbidity_meteotest(
        observer_latitude::Real,
        observer_longitude::Real,
        month_index::Int,
    )

[Bilinear interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation) of the Linke turbidity array, for the corrsponding latitude, longitude, and month index values. 

In the the Linke turbidity array, rows represent global latitudes from 90 to -90 degrees; columns represent global longitudes from -180 to 180; depth (third dimension) represents months of the year from January (1) to December (12). 

The array is based on the tabular data in [remund2009aerosol; Sec. 8](@cite).
"""
function linke_turbidity_meteotest(
    observer_latitude::Real,
    observer_longitude::Real,
    month_index::Int,
)
    i = size(lt_data, 1) * (latitude - LATITUDE_MIN) / LATITUDE_RANGE 
    j = size(lt_data, 2) * (longitude - LONGITUDE_MIN) / LONGITUDE_RANGE 

    i⁻, i⁺ = floor(i) |> Int, ceil(i) |> Int 
    j⁻, j⁺ = floor(j) |> Int, ceil(j) |> Int 

    # Index position between rounded row and column as [0, 1] fraction
    ν, τ = i - i⁻, j - j⁻

    # Bilinear interpolation of row-column data
    lt = ν * τ * LINKE_TURBIDITY_METEOTEST[i⁻, j⁻, month_index] + 
         ν * (1 - τ) * LINKE_TURBIDITY_METEOTEST[i⁻, j⁺, month_index] + 
         (1 - ν) * τ * LINKE_TURBIDITY_METEOTEST[i⁺, j⁻, month_index] + 
         (1 - ν) * (1 - τ) * LINKE_TURBIDITY_METEOTEST[i⁺, j⁺, month_index] 

    # The values within the Linke turbidity array are scaled by 20
    return lt / 20
end

"""
    linke_turbidity_kasten1996(
        absolute_airmass::Real, 
        precipitable_water::Real, 
        aod_bb::Real
    )

Calculate Linke turbidity using Kasten pyrheliometric formula.

Note that broadband aerosol optical depth (AOD) can be approximated by AOD measured at 700 nm according to Molineaux [molineaux1998equivalence](@cite). Bird and Hulstrom offer an alternate approximation using AOD measured at 380 nm and 500 nm.

Based on original implementation by Armel Oumbe. Note that these calculations are only valid for airmass less than 5 and precipitable water less than 5 cm.

# Arguments
- `absolute_airmass::Real` - Pressure-adjusted airmass.
- `precipitable_water::Real` - [cm] Precipitable water. 
- `aod_bb::Real` - broadband AOD. 
"""
function linke_turbidity_kasten1996(
    absolute_airmass::Real, 
    precipitable_water::Real, 
    aod_bb::Real
)
    # "From numerically integrated spectral simulations done with Modtran
    # (Berk, 1989), Molineaux (1998) obtained for the broadband optical depth
    # of a clean and dry atmospshere (fictitious atmosphere that comprises only
    # the effects of Rayleigh scattering and absorption by the atmosphere gases
    # other than the water vapor) the following expression"
    # - P. Ineichen (2008)
    δ_cda = -0.101 + 0.235absolute_airmass^(-0.16)
    # "and the broadband water vapor optical depth where pwat is the integrated
    # precipitable water vapor content of the atmosphere expressed in cm and am
    # the optical air mass. The precision of these fits is better than 1% when
    # compared with Modtran simulations in the range 1 < am < 5 and
    # 0 < pwat < 5 cm at sea level" - P. Ineichen (2008)
    δ_w = 0.112absolute_airmass^(-0.55) * precipitable_water^0.34
    # "Then using the Kasten pyrheliometric formula (1980, 1996), the Linke
    # turbidity at am = 2 can be written. The extension of the Linke turbidity
    # coefficient to other values of air mass was published by Ineichen and
    # Perez (2002)" - P. Ineichen (2008)
    return -(9.4 + 0.9absolute_airmass) * 
            log(exp(-absolute_airmass * (δ_cda + δ_w + aod_bb))) / absolute_airmass
end

"""
    pressure2altitude(pressure::Real)

Determine altitude, ``A`` [m], from site `pressure` ``P`` [Pa]. The following assumptions [psas2004altitudepressure](@cite) hold:

- Base pressure: 101325 Pa
- Temperature at zero altitude: 288.15 K
- Gravitational acceleration: 9.80665 m/s^2
- Lapse rate: -6.5E-3 K/m
- Gas constant for air: 287.053 J/(kg K)
- Relative Humidity: 0%

```math
    A = 44331.5 - 4946.62P^{0.190263}
```
"""
function pressure2altitude(pressure::Real)
    return 44331.5 - 4946.62pressure^0.190263
end

"""
    altitude2pressure(altitude::Real)

Determine site pressure, ``P`` [Pa], from altitude ``A`` [m]. It follows the same assumption of [`pressure2altitude`](@ref).

```math
    P = 100.0\\left(\\frac{44331.514 - A}{11880.516}\\right)^{\\frac{1}{0.1902632}}
```
"""
function altitude2pressure(altitude::Real)
    return 100.0((44331.514 - altitude) / 11880.516)^(1 / 0.1902632)
end

"""
    absolute_airmass(relative_airmass::Real, pressure::Real=101325.0)

Determine absolute (pressure-adjusted) airmass, ``a_\\text{abs}`` from relative airmass, ``a_\\text{reò}``, and pressure, ``P`` [Pa], [gueymard1993critical](@cite).

```math
    a_{\\text{abs}} = a_{\\text{rel}} \\frac{P}{P_{\\text{atm}}} 
```
where ``P_{\\text{atm}}`` is the atmospheric pressure.  
"""
function absolute_airmass(relative_airmass::Real, pressure::Real)
    return  relative_airmass * pressure / ATMOSPHERIC_PRESSURE
end

"""
    relative_airmass_simple(sun_apparent_elevation::Real)

Calculates the relative airmass as ``$VN_RELATIVE_AIRMASS = -\\sec $VN_TOPOCENTRIC_APPARENT_ELEVATION``, where $VN_TOPOCENTRIC_APPARENT_ELEVATION is the Sun's [topocentric_apparent_elevation](@ref) [deg]. Note that this returns `-Inf` at `sun_apparent_elevation=0`.
"""
function relative_airmass_simple(sun_apparent_elevation::Real)
    return -secd(sun_apparent_elevation)
end

"""
    relative_airmass_kasten1966(sun_apparent_elevation::Real)

Computes the relative airmass as proposed in [kasten1965new](@cite) as follows 
```math
$VN_RELATIVE_AIRMASS = \\frac{1}{-\\cos $VN_TOPOCENTRIC_APPARENT_ELEVATION + 0.15((3.885 + $VN_TOPOCENTRIC_APPARENT_ELEVATION)^(-1.253))}
```
where ``$VN_TOPOCENTRIC_APPARENT_ELEVATION`` is the Sun's [`topocentric_apparent_elevation`](@ref) [deg].
"""
function relative_airmass_kasten1966(sun_apparent_elevation::Real)
    e = sun_apparent_elevation
    return 1.0 / (-cosd(e) + 0.15((3.885 + e)^(-1.253)))
end

"""
    relative_airmass_youngirvine1967(sun_elevation::Real)

Computes the relative airmass as proposed in [young1967multicolor](@cite) as follows
```math
    $VN_RELATIVE_AIRMASS = -\\sec $VN_TOPOCENTRIC_ELEVATION (1 - 0.0012(\\sec $VN_TOPOCENTRIC_ELEVATION^2 + 1))
```
where ``$VN_TOPOCENTRIC_ELEVATION`` is the Sun's [`topocentric_elevation`](@ref) [deg].
"""
function relative_airmass_youngirvine1967(sun_elevation::Real)
    sec_elev = -secd(sun_elevation)
    return sec_elev * (1 - 0.0012(sec_elev^2 - 1))
end

"""
    relative_arimass_kastenyoung1989(sun_apparent_elevation::Real)

Computes the relative airmass as proposed in [kasten1989revised](@cite) as follows
```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{-\\cos $VN_TOPOCENTRIC_APPARENT_ELEVATION + 0.50572((6.07995 + $VN_TOPOCENTRIC_APPARENT_ELEVATION)^{-1.6364})}
```
where ``$VN_TOPOCENTRIC_APPARENT_ELEVATION`` is the Sun's [`topocentric_apparent_elevation`](@ref) [deg].
"""
function relative_arimass_kastenyoung1989(sun_apparent_elevation::Real)
    e = sun_apparent_elevation
    return 1.0 / (-cosd(e) + 0.50572((6.07995 + e)^(-1.6364)))
end

"""
    relative_airmass_gueymard1993(sun_apparent_elevation::Real)

Computes the relative airmass as proposed in [gueymard1993critical, gueymard1993development](@cite) as follows
```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{-\\cos $VN_TOPOCENTRIC_APPARENT_ELEVATION + 0.00176759(90 - $VN_TOPOCENTRIC_APPARENT_ELEVATION)((4.37515 + $VN_TOPOCENTRIC_APPARENT_ELEVATION)^{- 1.21563})}
```
where ``$VN_TOPOCENTRIC_APPARENT_ELEVATION`` is the Sun's [`topocentric_apparent_elevation`](@ref) [deg].
"""
function relative_airmass_gueymard1993(sun_apparent_elevation::Real)
    e = sun_apparent_elevation
    return 1.0 / (-cosd(e) + 0.00176759(90 - e)*((4.37515 + e)^(- 1.21563)))
end

"""
    relative_airmass_young1994(sun_elevation::Real)

Computes the relative airmass as proposed in [young1994air](@cite) as follows
```math
    $VN_RELATIVE_AIRMASS = \\frac{1.002432\\cos $VN_TOPOCENTRIC_ELEVATION^2 + 0.148386\\cos $VN_TOPOCENTRIC_ELEVATION + 0.0096467}{\\cos $VN_TOPOCENTRIC_ELEVATION^3 + 0.149864\\cos $VN_TOPOCENTRIC_ELEVATION^2 + 0.0102963\\cos $VN_RELATIVE_AIRMASS + 0.000303978}
```
where ``$VN_TOPOCENTRIC_ELEVATION`` is the Sun's [`topocentric_elevation`](@ref) [deg].
"""
function relative_airmass_young1994(sun_elevation::Real)
    cos_elev = -cosd(sun_elevation)
    return (1.002432cos_elev^2 + 0.148386cos_elev + 0.0096467) /
           (cos_elev^3 + 0.149864cos_elev^2 + 0.0102963cos_elev + 0.000303978)
end

"""
    relative_airmass_pickering2002(sun_apparent_elevation::Real)

Computes the relative airmass as proposed in [pickering2002southern](@cite) as follows
```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{\\sin($VN_TOPOCENTRIC_APPARENT_ELEVATION + 244.0 / (165 + 47.0$VN_TOPOCENTRIC_APPARENT_ELEVATION^{1.1}))}
```
where ``$VN_TOPOCENTRIC_ELEVATION`` is the Sun's [`topocentric_elevation`](@ref) [deg].
"""
function relative_airmass_pickering2002(sun_apparent_elevation::Real)
    e = sun_apparent_elevation
    return 1.0 / (sind(e + 244.0 / (165 + 47.0e^1.1)))
end

"""
    relative_airmass_gueymard2003(sun_apparent_elevation::Real)

Computes the relative airmass as proposed in [gueymard2003direct, gueymard2019clear](@cite) as follow as follows
```math
    $VN_RELATIVE_AIRMASS = -\\cos $VN_TOPOCENTRIC_APPARENT_ELEVATION + \\frac{0.48353(90 - $VN_TOPOCENTRIC_APPARENT_ELEVATION)^{0.095846}}{(6.741 + $VN_TOPOCENTRIC_APPARENT_ELEVATION)^{1.754}}
```
where ``$VN_TOPOCENTRIC_ELEVATION`` is the Sun's [`topocentric_elevation`](@ref) [deg].
"""
function relative_airmass_gueymard2003(sun_apparent_elevation::Real)
    e = sun_apparent_elevation
    return -cosd(e) + 0.48353(90 - e)^0.095846 / (6.741 + e)^1.754
end

"""
    precipitable_water_gueymard94(
        air_temperature::Real, 
        relative_humidity::Real
    )

Calculates precipitable water (cm) from ambient air temperature [C] and relatively humidity [%] using an empirical model. The accuracy of this method is approximately 20% for moderate PW (1-3 cm) and less accurate otherwise.

The model was developed by expanding Eq. 1 in [gueymard1994analysis; Eq. 1](@cite) and using [gueymard1994analysis; Eq. 2](@cite):

```math
    \\begin{align*}
        Pw = 0.1 H_v \\rho_v \\\\
        \\rho_v = 216.7 R_H e_s
    \\end{align*}
```

``Pw`` is the precipitable water (cm), ``H_v`` is the apparent water vapor scale height [km] and :math:``\\rho_v`` is the surface water vapor density [g/m^3]. The expression for ``H_v`` results from [gueymard1994analysis; Eq. 4](@cite):

```math
    H_v = 0.4976 + 1.5265 \\frac{T}{273.15}
        + \\exp \\left(13.6897 \\frac{T}{273.15}
        - 14.9188 \\left( \\frac{T}{273.15} \\right)^3 \\right)
```

In the expression for ``\\rho_v``, ``e_s`` is the saturation water vapor pressure [millibar]. The expression for ``e_s`` results from [keogh2004accurate; Eq. 1](@cite)

```math
   e_s = \\exp  \\left(22.330 - 49.140 \\frac{100}{T} -
       10.922 \\left(\\frac{100}{T}\\right)^2 -
       0.39015 \\frac{T}{100} \\right)
```
"""
function precipitable_water_gueymard94(
    air_temperature::Real, 
    relative_humidity::Real
)
    T = air_temperature + 273.15  # Convert to Kelvin
    rh = relative_humidity

    θ = T / 273.15

    # Eq. 1 from Keogh and Blakers
    pw = 0.1(0.4976 + 1.5265θ + exp(13.6897θ - 14.9188θ^3)) * 
         (216.7rh/(100T)*exp(22.330 - 49.140(100/T) - 10.922(100/T)^2 - 0.39015T/100))

    return maximum(pw, 0.1)
end

"""
    dew_temperature2relative_humidity_oke2018(
        air_temperature::Real,
        dew_temperature::Real
    )

Calculate relative humidity [%] from dewpoint temperature [°C] using the Magnus equation. Air temperature is expressed in °C. The Magnus equation coefficients (A, B, C) default to the values recommended by [oke2018guide](@cite).
"""
function dew_temperature2relative_humidity_oke2018(
    air_temperature::Real,
    dew_temperature::Real
)

    # Calculate vapor pressure (e) and saturation vapor pressure (es)
    e = 6.112exp(17.62temp_air / (243.12 + temp_air))
    es = 6.112exp(17.62temp_dew / (243.12 + temp_dew))

    # Calculate relative humidity as percentage
    return 100es / e
end

"""
    relative_humidity2dew_temperature_oke2018(
        air_temperature::Real,
        relative_humidity::Real
    )

Calculate dewpoint temperature using the Magnus equation. It corresponds to the inverse of [`dew_temperature2relative_humidity_oke2018`](@ref).
"""
function relative_humidity2dew_temperature_oke2018(
    air_temperature::Real, 
    relative_humidity::Real
)
    # Calculate the term inside the log
    # From RH = 100 * (es/e), we get es = (RH/100) * e
    # Substituting the Magnus equation and solving for dewpoint

    # First calculate ln(es/A)
    ln_term = 6.112temp_air / (coeff[2] + temp_air) + log(relative_humidity/100)
    
    # Then solve for dewpoint
    return coeff[2] * ln_term / (coeff[1] - ln_term)
end
