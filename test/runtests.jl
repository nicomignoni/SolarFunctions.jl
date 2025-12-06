using Test, Dates, SolarFunctions

print_error(name, test, correct) = 
    println(
        "$(name):
        - test: $(test)
        - correct: $(correct)
        - error: $(100(test - correct)/correct |> abs) %\n"
    )

@testset "solar position" begin
    # Test data
    latitude = 39.742476 
    longitude = -105.1786
    altitude = 1830.14
    temperature = 11.0
    pressure = 820
    tz = -7
    time = DateTime(2003, 10, 17, 12 - tz, 30, 30) 
    ΔT = 67


    correct = (
        jd = 2452930.312847,
        jc = 0.037928, 
        L = 24.0182616917,
        B = -0.0001011219,
        R = 0.9965422974,
        Θ = 204.0182616917,
        Δψ = -0.00399840,
        ϵ = 23.440465,
        α = 202.22741,
        H = 11.105900,
        ᾱ = 202.22704,
        β = 0.0001011219,
        Δϵ = 0.00166657,
        λ = 204.0085519281,
        δ = -9.31434,
        H̄ = 11.10629,
        δ̄ = -9.316179,
        Φ = 194.34024,
        η = 39.88838
    )

    jd = SolarFunctions.datetime2julian(time)
    jed = SolarFunctions.julian_ephemeris_day(jd, ΔT)
    jc = SolarFunctions.julian_century(jd)
    jec = SolarFunctions.julian_ephemeris_century(jed)
    jem = SolarFunctions.julian_ephemeris_millenium(jec)

    print_error("Julian day", jd, correct.jd)
    print_error("Julian century", jc, correct.jc)

    L = SolarFunctions.heliocentric_longitude(jem)
    B = SolarFunctions.heliocentric_latitude(jem)
    R = SolarFunctions.heliocentric_radius(jem)
    Θ = SolarFunctions.geocentric_longitude(L)
    β = SolarFunctions.geocentric_latitude(B)
    
    print_error("Helioc. longitude", L, correct.L)
    print_error("Helioc. latitude", B, correct.B)
    print_error("Helioc. radius", R, correct.R)
    print_error("Geoc. longitude", Θ, correct.Θ)
    print_error("Geoc. longitude", β, correct.β)

    W = SolarFunctions.nutation_coefficients(jec)
    Δψ = SolarFunctions.nutation_longitude(jec, W)
    Δϵ = SolarFunctions.nutation_obliquity(jec, W)
    ϵ₀ = SolarFunctions.mean_elliptic_obliquity(jem)
    ϵ = SolarFunctions.elliptic_obliquity(ϵ₀, Δϵ)

    print_error("Nutat. longitude", Δψ, correct.Δψ)
    print_error("Nutat. obliquity", Δϵ, correct.Δϵ)
    print_error("Ellip. obliquity", ϵ, correct.ϵ)

    Δτ = SolarFunctions.aberration_correction(R)
    λ = SolarFunctions.apparent_sun_longitude(Θ, Δψ, Δτ)
    ν₀ = SolarFunctions.mean_sidereal_greenwich_time(jd, jc)
    ν = SolarFunctions.apparent_sidereal_greenwich_time(ν₀, Δψ, ϵ)
    α = SolarFunctions.geocentric_sun_ascension(λ, ϵ, β)
    δ = SolarFunctions.geocentric_sun_declination(λ, ϵ, β)
    H = SolarFunctions.observer_local_hour(longitude, ν, α)
    ξ = SolarFunctions.sun_equatorial_horizontal_parallax(R)

    print_error("App. Sun longitude", λ, correct.λ)
    print_error("Geoc. Sun ascension", α, correct.α)
    print_error("Geoc. Sun declination", δ, correct.δ)
    print_error("Obser. hour angle", H, correct.H)

    u = SolarFunctions.reduced_observer_latitude(latitude)
    x = SolarFunctions.radial_distance_equatorial_plane(latitude, altitude, u) 
    y = SolarFunctions.radial_distance_rotational_axis(latitude, altitude, u)
    Δα = SolarFunctions.sun_ascension_parallax(x, ξ, δ, H)
    δ̄ = SolarFunctions.topocentric_sun_declination(x, y, ξ, δ, H, Δα)
    ᾱ = SolarFunctions.topocentric_sun_ascension(α, Δα)
    H̄ = SolarFunctions.topocentric_local_hour(H, Δα)
    η₀ = SolarFunctions.topocentric_apparent_elevation(latitude, δ̄, H̄)
    η̂₀ = SolarFunctions.topocentric_elevation_correction(temperature, pressure, η₀)
    η = SolarFunctions.topocentric_elevation(η₀, η̂₀)
    Γ = SolarFunctions.topocentric_astronomical_azimuth(latitude, δ̄, H̄)
    Φ = SolarFunctions.topocentric_azimuth(Γ)

    print_error("Topoc. Sun declination", δ̄, correct.δ̄)
    print_error("Topoc. Sun ascension", ᾱ, correct.ᾱ)
    print_error("Topoc. hour angle", H̄, correct.H̄)
    print_error("Topoc. elevation", η, correct.η)
end

@testset "clearsky" begin
    sun_apparent_elevation = 30.4 
    observer_altitude = 0.0
    absolute_airmass = 1.0
    linke_turbidity = 2.1
    extraterrestial_radiation = 1364.0
    perez_enanchement = false

    correct = (
        ineichen = (
            dni = 1007.606890209482,
            dhi = 42.47213240572813,
            ghi = 552.3552398128525,
        ), 
        haurwitz = (
            ghi = 494.477044,
        )
    )

    # Ineichen
    dni, dhi, ghi = SolarFunctions.clearsky_ineichen(
        sun_apparent_elevation,
        observer_altitude,
        absolute_airmass,
        linke_turbidity,
        extraterrestial_radiation,
        perez_enanchement
    )
    print_error("Ineichen (DNI)", dni, correct.ineichen.dni)
    print_error("Ineichen (DHI)", dhi, correct.ineichen.dhi)
    print_error("Ineichen (GHI)", ghi, correct.ineichen.ghi)

    # Haurwitz
    ghi = SolarFunctions.clearsky_haurwitz(sun_apparent_elevation)
    print_error("Haurwitz (GHI)", ghi, correct.haurwitz.ghi)
end
