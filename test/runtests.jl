using Test, Dates, Serialization, Irradia

@testset "solar position" begin
    # Test data
    latitude = 39.742476 
    longitude = -105.1786
    altitude = 1830.14
    temperature = 11.0
    pressure = 820
    tz = -7
    time = DateTime(2003, 10, 17, 12 - tz, 30, 30) 

    correct = Dict(
        "jd" => 2452930.312847,
        "jc" => 0.037928, 
        "L" => 24.0182616917,
        "B" => -0.0001011219,
        "R" => 0.9965422974,
        "θ" => 204.0182616917,
        "Δψ" => -0.00399840,
        "ϵ" => 23.440465,
        "α" => 202.22741,
        "H" => 11.105900,
        "ᾱ" => 202.22704,
        "β" => 0.0001011219,
        "Δϵ" => 0.00166657,
        "λ" => 204.0085519281,
        "δ" => -9.31434,
        "H̄" => 11.10629,
        "δ̄" => -9.316179,
        "Φ" => 194.34024,
        "e" => 39.88838
    )

    test = Dict()
    test["jd"] = datetime2julian(time)
    test["jed"] = julian_ephemeris_day(test["jd"])
    test["jc"] = julian_century(test["jd"])
    test["jec"] = julian_ephemeris_century(test["jed"])
    test["jem"] = julian_ephemeris_millenium(test["jec"])
    test["L"] = heliocentric_polynomial(Irradia.EARTH_PERIODIC_TERMS.L, test["jem"]) |> rad2deg |> Irradia.mod360
    test["B"] = heliocentric_polynomial(Irradia.EARTH_PERIODIC_TERMS.B, test["jem"]) |> rad2deg
    test["R"] = heliocentric_polynomial(Irradia.EARTH_PERIODIC_TERMS.R, test["jem"])
    test["θ"] = Irradia.mod360(test["L"] + 180)
    test["β"] = -test["B"]
    test["X"] = [
        mean_moon_elongation_from_sun(test["jec"]),
        mean_sun_anomaly(test["jec"]),
        mean_moon_anomaly(test["jec"]),
        moon_latitude_argument(test["jec"]),
        moon_ascdending_longitude(test["jec"])
    ]
    test["Z"] = Irradia.EARTH_PERIODIC_TERMS.Y * test["X"] # (in degrees)
    test["Δψ"] = nutation_longitude(test["jec"], test["Z"], Irradia.EARTH_PERIODIC_TERMS.ψ)
    test["Δϵ"] = nutation_obliquity(test["jec"], test["Z"], Irradia.EARTH_PERIODIC_TERMS.ϵ)
    test["ϵ"] = elliptic_obliquity(test["jem"], test["Δϵ"])
    test["Δτ"] = aberration_correction(test["R"])
    test["λ"] = apparent_sun_longitude(test["θ"], test["Δψ"], test["Δτ"])
    test["ν"] = apparent_sidereal_greenwich_time(test["jd"], test["jc"], test["Δψ"], test["ϵ"])
    test["α"] = geodetic_sun_ascension(test["λ"], test["ϵ"], test["β"])
    test["δ"] = geodetic_sun_declination(test["λ"], test["ϵ"], test["β"])
    test["H"] = observer_local_hour_angle(test["ν"], longitude, test["α"])
    test["ξ"] = equatorial_horizontal_sun_parallax(test["R"])
    test["u"] = atand(0.99664719 * tand(latitude)) # (in degrees)
    test["x"] = cosd(test["u"]) + altitude * cosd(latitude) / 6378140.0
    test["y"] = 0.99664719 * sind(test["u"]) + altitude * sind(latitude) / 6378140.0
    test["Δα"] = sun_ascension_parallax(test["x"], test["y"], test["ξ"], test["δ"], test["H"])
    test["δ̄"] = topocentric_sun_declination(test["x"], test["y"], test["ξ"], test["δ"], test["H"], test["Δα"])
    test["ᾱ"] = topocentric_sun_ascension(test["α"], test["Δα"])
    test["H̄"] = topocentric_local_hour_angle(test["H"], test["Δα"])
    test["e"] = topocentric_elevation_angle(latitude, temperature, pressure, test["δ̄"], test["H̄"])
    test["Γ"] = topocentric_astronomical_azimuth_angle(latitude, test["δ̄"], test["H̄"])
    test["Φ"] = topocentric_azimuth_angle(test["Γ"])

    for key in keys(correct)
        error = (correct[key] - test[key]) / test[key] |> abs
        println("$key: $(round(100error, digits=5)) %")
    end
end
