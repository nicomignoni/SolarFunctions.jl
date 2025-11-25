using Dates, Serialization

# The .h5 file 'LinkeTurbidities.h5' contains a single 2160 x 4320 x 12
# matrix of type uint8 called 'LinkeTurbidity'.  To determine the Linke
# turbidity for a position on the Earth's surface for a given month do the
# following: LT = LinkeTurbidity(LatitudeIndex, LongitudeIndex, month).

# The nodes of the grid are 5' (1/12=0.0833[arcdeg]) apart.
# From Section 8 of Aerosol optical depth and Linke turbidity climatology
# http://www.meteonorm.com/images/uploads/downloads/ieashc36_report_TL_AOD_climatologies.pdf
# 1st row: 89.9583 S, 2nd row: 89.875 S
# 1st column: 179.9583 W, 2nd column: 179.875 W
const LINKE_TURBIDITY = open("data/linke-turbidity.jld", "r") do io;
    deserialize(io)
end

# Bounds (in degrees) for latitude and longitude
const LATITUDE_MAX = 90
const LATITUDE_MIN = -90
const LONGITUDE_MAX = 180
const LONGITUDE_MIN = -180 

const latitude_range = LATITUDE_MAX - LATITUDE_MIN
const longitude_range = LONGITUDE_MAX - LONGITUDE_MIN

"""
    linke_turbidity(lat::Real, lon::Real, month::Int, lt_data::Array)

Interpolates the Like turbidity for the corrsponding latitude, longitude, and month values.

# Arguments
- lat::Real - [degrees] the observer geographical latitude
- lon::Real - [degrees] the observer geographical longitude
- month::Int - Month expressed in a 1-12 range (i.e., 1 for January, 5 for May)
- `lt_data::Array` - Linke turbidity array. Rows represent global latitudes from 90 
to -90 degrees; columns represent global longitudes from -180 to 180; depth 
(third dimension) represents months of the year from January (1) to December (12). 
"""
function linke_turbidity(
    latitude::Real,
    longitude::Real,
    month::Int,
    lt_data::Array
)
    i = size(turbidity_data, 1) * (latitude - LATITUDE_MIN) / latitude_range 
    j = size(turbidity_data, 2) * (longitude - LONGITUDE_MIN) / longitude_range 

    i⁻, i⁺ = floor(i) |> Int, ceil(i) |> Int 
    j⁻, j⁺ = floor(j) |> Int, ceil(j) |> Int 

    # Index position between rounded row and column as [0, 1] fraction
    ν, τ = i - i⁻, j - j⁻

    # Bilinear interpolation of row-column data
    lt = ν * τ * lt_data[i⁻, j⁻, month] + 
         ν * (1 - τ) * lt_data[i⁻, j⁺, month] + 
         (1 - ν) * τ * lt_data[i⁺, j⁻, month] + 
         (1 - ν) * (1 - τ) * lt_data[i⁺, j⁺, month] 

    # The numbers within the matrix are 20 * Linke Turbidity
    return lt / 20
end

"""
    ineichen(apparent_zenith, airmass_absolute, linke_turbidity, altitude=0, dni_extra=1364.0, perez_enhancement::Bool=false)

Determine clear sky GHI, DNI, and DHI from Ineichen/Perez model.

Implements the Ineichen and Perez clear sky model for global
horizontal irradiance (GHI), direct normal irradiance (DNI), and
calculates the clear-sky diffuse horizontal (DHI) component as the
difference between GHI and DNI*cos(zenith) as presented in [^1] [^2]. A
report on clear sky models found the Ineichen/Perez model to have
excellent performance with a minimal input data set [^3].

Default values for monthly Linke turbidity provided by SoDa [^4], [^5].

# Arguments
- `apparent_zenith::Real` - Refraction corrected solar zenith angle in degrees.
- `airmass_absolute::Real` - Pressure corrected airmass.
- `linke_turbidity::Real` - Linke Turbidity.
- `altitude::Real = 0` - Altitude above sea level in meters.
- `dni_extra::Real=1364` - Extraterrestrial irradiance. The units of ``dni_extra`` determine the units of the output.
- `perez_enhancement::Bool=false` - Controls if the Perez enhancement factor should be applied. Setting to `true` [may 
produce spurious results](https://github.com/pvlib/pvlib-python/issues/435) for times when the Sun is near the horizon 
and the airmass is high.

# Returns
- `dni::Real` - direct normal irradiance 
- `dhi::Real` - direct horizontal irradiance 
- `ghi::Real` - global horizoantal irradiance

[^1]: P. Ineichen and R. Perez, "A New airmass independent formulation for
the Linke turbidity coefficient", Solar Energy, vol 73, pp. 151-157,
2002.

[^2]: R. Perez et. al., "A New Operational Model for Satellite-Derived
Irradiances: Description and Validation", Solar Energy, vol 73, pp.
307-317, 2002.

[^3]: M. Reno, C. Hansen, and J. Stein, "Global Horizontal Irradiance
Clear Sky Models: Implementation and Analysis", Sandia National
Laboratories, SAND2012-2389, 2012.

[^4]: https://www.soda-pro.com/help/general-knowledge/linke-turbidity-factor
(accessed February 2, 2024).

[^5]: J. Remund, et. al., "Worldwide Linke Turbidity Information", Proc.
ISES Solar World Congress, June 2003. Goteborg, Sweden.
"""
function ineichen(
    apparent_zenith,
    airmass_absolute,
    linke_turbidity,
    altitude=0,
    dni_extra=1364.0,
    perez_enhancement::Bool=false
)
    cos_zenith = maximum(cosd(apparent_zenith), 0)

    tl = linke_turbidity

    fh1 = exp(-altitude / 8000.0)
    fh2 = exp(-altitude / 1250.0)
    cg1 = 5.09e-05 * altitude + 0.868
    cg2 = 3.92e-05 * altitude + 0.0387

    ghi = exp(-cg2 * airmass_absolute * (fh1 + fh2 * (tl - 1)))

    # https://github.com/pvlib/pvlib-python/issues/435
    if perez_enhancement
        ghi *= exp(0.01 * airmass_absolute^1.8)
    end

    # use maximum to map airmass nans to 0s. multiply and divide by tl to
    # reinsert tl nans
    ghi = cg1 * dni_extra * cos_zenith * tl / tl * maximum(ghi, 0)

    # From [1] (Following [2] leads to 0.664 + 0.16268 / fh1)
    # See https://github.com/pvlib/pvlib-python/pull/808
    b = 0.664 + 0.163 / fh1
    # BncI = "normal beam clear sky radiation"
    bnci = b * exp(-0.09 * airmass_absolute * (tl - 1))
    bnci = dni_extra * maximum(bnci, 0)

    # "empirical correction" SE 73, 157 & SE 73, 312.
    bnci_2 = ((1 - (0.1 - 0.2 * exp(-tl)) / (0.1 + 0.882 / fh1)) / cos_zenith)
    bnci_2 = ghi * minimum(maximum(bnci_2, 0), 1e20)

    dni = minimum(bnci, bnci_2)

    dhi = ghi - dni * cos_zenith

    return dni, dhi, ghi
end
